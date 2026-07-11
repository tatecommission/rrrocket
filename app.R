library(shiny)
library(bslib)
library(plotly)
library(ggplot2)
library(leaflet)
library(leaflet.extras)
library(jsonlite)
library(xml2)

# =============================================================================
# RRRocket 3D
#
# Physics modeled on:
#   Niskanen, S. (2013). OpenRocket technical documentation, v13.05.
#   Barrowman, J. (1967). The Practical Calculation of the Aerodynamic
#     Characteristics of Slender Finned Vehicles.
#   Hoerner, S. (1965). Fluid-Dynamic Drag.
#
# Equation numbers below (eq. 3.xx / 4.xx / B.x) refer to Niskanen 2013.
#
# Integration: Runge-Kutta 4 (eq. 4.20-4.21).
# Rotational model: pitch-plane rigid body. State carries the body axis unit
#   vector and the angular velocity vector, so angle of attack, weathercocking,
#   gravity turn and divergence of unstable rockets all EMERGE from the
#   normal force acting at the CP rather than being applied as a heuristic.
# =============================================================================

unit_choices_length <- c("mm", "cm", "in", "ft", "m")
unit_choices_mass   <- c("g", "oz", "kg", "lb")
unit_choices_speed  <- c("m/s", "mph", "km/h", "knots")

g0      <- 9.80665
R_air   <- 287.058
gamma_a <- 1.4
Cd_chute <- 0.80          # Hoerner parachute Cd (Niskanen sec. 4.2.5)
m_to_ft  <- 3.28084
N_to_lbf <- 0.224809

# Safety thresholds
V_DEPLOY_WARN <- 20       # m/s -- deployment above this risks a zippered tube
V_RAIL_WARN   <- 15       # m/s -- minimum safe rail-exit speed

# ---- unit conversions -------------------------------------------------------
to_meters <- function(val, unit) {
  if (!isTruthy(val)) return(0)
  switch(unit, "mm"=val/1000, "cm"=val/100, "in"=val*0.0254,
         "ft"=val*0.3048, "m"=val, val/1000)
}
to_kg <- function(val, unit) {
  if (!isTruthy(val)) return(0)
  switch(unit, "g"=val/1000, "oz"=val*0.0283495,
         "kg"=val, "lb"=val*0.453592, val/1000)
}
to_ms <- function(val, unit) {
  if (!isTruthy(val)) return(0)
  switch(unit, "m/s"=val, "mph"=val*0.44704,
         "km/h"=val/3.6, "knots"=val*0.514444, val)
}
from_meters <- function(si, unit) {
  switch(unit, "mm"=si*1000, "cm"=si*100, "in"=si/0.0254,
         "ft"=si/0.3048, "m"=si, si*1000)
}
from_kg <- function(si, unit) {
  switch(unit, "g"=si*1000, "oz"=si/0.0283495,
         "kg"=si, "lb"=si/0.453592, si*1000)
}
from_ms <- function(si, unit) {
  switch(unit, "m/s"=si, "mph"=si/0.44704,
         "km/h"=si*3.6, "knots"=si/0.514444, si)
}

# ---- small vector helpers ---------------------------------------------------
cross3 <- function(a, b) c(a[2]*b[3]-a[3]*b[2],
                           a[3]*b[1]-a[1]*b[3],
                           a[1]*b[2]-a[2]*b[1])
vnorm  <- function(a) sqrt(sum(a*a))
unitv  <- function(a) { n <- vnorm(a); if (n < 1e-12) c(0,0,1) else a/n }

# ---- ISA troposphere (Niskanen sec. 4.1.1) ----------------------------------
isa_atm <- function(z) {
  z   <- max(z, 0)
  T   <- 288.15 - 0.0065 * z
  rho <- 1.225 * (T / 288.15)^4.2561
  list(rho = rho, a = sqrt(gamma_a * R_air * T))
}

# ---- Skin friction (eqs. 3.78-3.84) -----------------------------------------
# Fully turbulent boundary layer assumed (Niskanen sec. 3.4.1).
skin_friction_cf <- function(velocity, char_length, mach, Rs = 60e-6) {
  nu    <- 1.461e-5
  Re    <- max(velocity * char_length / nu, 1)
  Rcrit <- 51 * (Rs / char_length)^(-1.039)      # eq. 3.79
  
  Cf <- if (Re < 1e4) {
    1.48e-2                                       # eq. 3.81 low-Re floor
  } else if (Re < Rcrit) {
    1 / (1.50 * log(Re) - 5.6)^2                  # eq. 3.78
  } else {
    0.032 * (Rs / char_length)^0.2                # eq. 3.80
  }
  
  if (mach < 1.0) {
    Cf_c <- Cf * (1.0 - 0.1 * mach^2)             # eq. 3.82
  } else {
    Cf_c <- Cf / (1.0 + 0.15 * mach^2)^0.58       # eq. 3.83
    if (Re >= Rcrit) {
      # eq. 3.84, never below the turbulent value
      Cf_c <- max(Cf_c, Cf / (1.0 + 0.18 * mach^2))
    }
  }
  Cf_c
}

# ---- Nose pressure drag Mach scaling ----------------------------------------
mach_cd_factor <- function(M) {
  if      (M < 0.8) 1 / sqrt(max(1 - M^2, 0.01))
  else if (M < 1.0) 1.6667 + 3.6667 * (M - 0.8)
  else if (M < 2.0) 2.4 - 0.8 * (M - 1.0)
  else              1.6
}

# ---- Base drag (eq. 3.94) ---------------------------------------------------
live_base_drag <- function(M) {
  if (M < 1.0) 0.12 + 0.13 * M^2
  else         0.25 / M
}

# ---- Stagnation pressure coefficient (eqs. B.1, B.2) ------------------------
# Used for launch-lug / rail-pin parasitic drag.
cd_stag <- function(M) {
  qratio <- if (M < 1) {
    1 + M^2/4 + M^4/40
  } else {
    1.84 - 0.76/M^2 + 0.166/M^4 + 0.035/M^6
  }
  0.85 * qratio
}

# ---- Angle-of-attack axial drag scaling (Niskanen sec. 3.4.7) ---------------
# 1.0 at 0 deg, 1.3 at 17 deg, 0 at 90 deg; zero derivative at all three.
# Implemented as two smoothstep segments, which satisfies exactly those
# constraints and is monotone (no polynomial overshoot).
aoa_drag_factor <- function(alpha_rad) {
  a <- abs(alpha_rad)
  a <- min(a, pi - a)                 # symmetric: tail-first behaves like nose-first
  deg <- a * 180 / pi
  smooth <- function(t) { t <- min(max(t, 0), 1); 3*t^2 - 2*t^3 }
  if (deg <= 17) {
    1 + 0.3 * smooth(deg / 17)
  } else {
    1.3 * (1 - smooth((deg - 17) / (90 - 17)))
  }
}

# ---- Fin-fin interference (eq. 3.54) ----------------------------------------
fin_interference <- function(N) {
  if      (N <= 4) 1.000
  else if (N == 5) 0.948
  else if (N == 6) 0.913
  else if (N == 7) 0.862      # interpolated between N=6 and N=8
  else if (N == 8) 0.810
  else             0.750
}

# ---- Transonic fin CP position (eqs. 3.35-3.36) -----------------------------
# Supersonic limit, M > 2:  Xf/cbar = (A*beta - 0.67) / (2*A*beta - 1)
fin_cp_super  <- function(M, A) { b <- sqrt(M^2 - 1); (A*b - 0.67) / (2*A*b - 1) }
fin_cp_super_d <- function(M, A) {
  b <- sqrt(M^2 - 1)
  0.34 * A / (2*A*b - 1)^2 * (M / b)
}
# Quintic p(M) on [0.5, 2] with the six constraints of eq. 3.36:
#   p(0.5)=0.25, p'(0.5)=0, p(2)=f(2), p'(2)=f'(2), p''(2)=0, p'''(2)=0
fin_cp_poly <- function(A) {
  M0 <- 0.5; M1 <- 2
  rv  <- function(M) M^(0:5)
  rd1 <- function(M) c(0, 1, 2*M, 3*M^2, 4*M^3, 5*M^4)
  rd2 <- function(M) c(0, 0, 2, 6*M, 12*M^2, 20*M^3)
  rd3 <- function(M) c(0, 0, 0, 6, 24*M, 60*M^2)
  Amat <- rbind(rv(M0), rd1(M0), rv(M1), rd1(M1), rd2(M1), rd3(M1))
  bvec <- c(0.25, 0, fin_cp_super(M1, A), fin_cp_super_d(M1, A), 0, 0)
  tryCatch(as.numeric(solve(Amat, bvec)),
           error = function(e) c(0.25, 0, 0, 0, 0, 0))
}
fin_cp_frac <- function(M, A, coefs) {
  fr <- if (M <= 0.5)      0.25
  else if (M >= 2.0) fin_cp_super(M, A)
  else               sum(coefs * M^(0:5))
  min(max(fr, 0.05), 1.0)   # keep the CP physically on the chord
}

# =============================================================================
# GEOMETRY / STATIC AERODYNAMICS
# Everything Mach-independent is precomputed once here. Mach-dependent
# quantities (CNa, CP, CD0) are recomputed every RK4 stage by aero_coeffs().
# =============================================================================
compute_aero <- function(nose_type, nose_length, body_length,
                         bt_diameter, fin_count, fin_root,
                         fin_tip, fin_span, fin_sweep, fin_pos,
                         cg_measured,
                         lug_length = 0, lug_od = 0, lug_id = 0,
                         ref_velocity = 50) {
  
  bt_radius <- bt_diameter / 2
  Aref      <- pi * bt_radius^2
  N         <- fin_count
  
  # --- nose (Barrowman) ---
  Xcp_nose <- switch(nose_type,
                     conical   = (2/3)  * nose_length,
                     ogive     = 0.466  * nose_length,
                     parabolic = 0.5    * nose_length)
  CNa_nose <- 2.0
  
  # --- fin planform ---
  ct_sum   <- fin_root + fin_tip
  Afin_one <- 0.5 * ct_sum * fin_span                       # one side, one fin
  A_aspect <- 2 * fin_span^2 / Afin_one                     # eq. 3.35 aspect ratio
  # Mean aerodynamic chord (eq. 3.30 evaluated for a trapezoid)
  cbar     <- (2/3) * (fin_root^2 + fin_tip^2 + fin_root*fin_tip) / ct_sum
  # Leading-edge position of the MAC (eq. 3.32 -> Barrowman's Xt term)
  xmac_le  <- (fin_sweep / 3) * (fin_root + 2*fin_tip) / ct_sum
  # Midchord sweep angle, needed by eq. 3.40
  dx_mid   <- fin_sweep + fin_tip/2 - fin_root/2
  cosGc    <- fin_span / sqrt(fin_span^2 + dx_mid^2)
  
  Kfb      <- 1 + bt_radius / (fin_span + bt_radius)        # eq. 3.56
  fin_int  <- fin_interference(N)                           # eq. 3.54
  cp_coefs <- fin_cp_poly(A_aspect)
  
  # Static (subsonic) fin CP, also used as the fin lever arm for pitch damping
  Xcp_fin0 <- fin_pos + xmac_le + 0.25 * cbar               # eq. 3.34
  
  # --- nose pressure drag ---
  ha <- atan(bt_radius / nose_length)
  Cd_nose_pressure <- switch(nose_type,
                             conical   = 0.8 * sin(ha)^2,
                             ogive     = 0.5 * sin(ha)^2,
                             parabolic = 0.3 * (bt_diameter / nose_length)^2)
  
  rocket_length <- nose_length + body_length
  fB            <- rocket_length / bt_diameter
  nose_slant    <- sqrt(nose_length^2 + bt_radius^2)
  Awet_nose     <- pi * bt_radius * nose_slant
  Awet_body     <- pi * bt_diameter * body_length
  Awet_fins     <- 2 * N * Afin_one
  
  # --- launch lug parasitic drag (eqs. 3.95, 3.96) ---
  if (isTruthy(lug_length) && lug_length > 0 && lug_od > 0) {
    r_ext <- lug_od / 2
    r_int <- min(lug_id, lug_od) / 2
    ld    <- lug_length / lug_od
    lug_k <- max(1.3 - 0.3 * ld, 1)                                   # eq. 3.95
    A_lug <- pi*r_ext^2 - pi*r_int^2 * max(1 - ld, 0)                 # eq. 3.96
  } else {
    lug_k <- 0; A_lug <- 0
  }
  
  aero <- list(
    nose_type = nose_type, nose_length = nose_length, body_length = body_length,
    bt_diameter = bt_diameter, bt_radius = bt_radius, Aref = Aref,
    fin_count = N, fin_root = fin_root, fin_tip = fin_tip,
    fin_span = fin_span, fin_sweep = fin_sweep, fin_pos = fin_pos,
    Afin_one = Afin_one, A_aspect = A_aspect, cbar = cbar,
    xmac_le = xmac_le, cosGc = cosGc, Kfb = Kfb, fin_int = fin_int,
    cp_coefs = cp_coefs, Xcp_fin0 = Xcp_fin0,
    CNa_nose = CNa_nose, Xcp_nose = Xcp_nose,
    Cd_nose_pressure = Cd_nose_pressure,
    Awet_nose = Awet_nose, Awet_body = Awet_body, Awet_fins = Awet_fins,
    fB = fB, rocket_length = rocket_length,
    lug_k = lug_k, A_lug = A_lug,
    cd_scale = 1.0
  )
  
  # ---- display-only reference values (M = 0, ref_velocity) ----
  ref <- aero_coeffs(aero, ref_velocity, 0)
  Cf_ref     <- skin_friction_cf(ref_velocity, rocket_length, 0)
  Cf_fin_ref <- skin_friction_cf(ref_velocity, cbar,          0)
  
  aero$CNa_total <- ref$CNa
  aero$CP        <- ref$CP
  aero$Cd_nose   <- Cd_nose_pressure + Cf_ref * Awet_nose / Aref
  aero$Cd_body   <- Cf_ref * (1 + 2/fB) * Awet_body / Aref
  aero$Cd_fins   <- Cf_fin_ref * Awet_fins / Aref
  aero$Cd_base   <- 0.12
  aero$Cd_lug    <- lug_k * cd_stag(0) * A_lug / Aref
  aero$Cd        <- ref$CD0                 # nominal CD0 at 50 m/s, M = 0
  aero$stability_margin <- (ref$CP - cg_measured) / bt_diameter
  aero
}

# ---- Mach-/velocity-dependent coefficients (called every RK4 stage) ---------
aero_coeffs <- function(aero, vel, M) {
  # --- fin normal force, eq. 3.40 (reduces to Barrowman at M = 0) ---
  beta  <- max(sqrt(abs(1 - M^2)), 1e-3)
  denom <- 1 + sqrt(1 + (beta * aero$fin_span^2 /
                           (aero$Afin_one * aero$cosGc))^2)
  CNa1    <- 2*pi * (aero$fin_span^2 / aero$Aref) / denom
  CNa_fin <- aero$Kfb * aero$fin_int * (aero$fin_count / 2) * CNa1
  
  # --- fin CP marches aft transonically (eqs. 3.35-3.36) ---
  frac    <- fin_cp_frac(M, aero$A_aspect, aero$cp_coefs)
  Xcp_fin <- aero$fin_pos + aero$xmac_le + frac * aero$cbar
  
  CNa <- aero$CNa_nose + CNa_fin
  CP  <- (aero$CNa_nose * aero$Xcp_nose + CNa_fin * Xcp_fin) / CNa
  
  # --- zero-AoA axial drag, CD0 (eq. 3.97) ---
  vs      <- max(vel, 0.1)
  Cd_pres <- aero$Cd_nose_pressure * mach_cd_factor(M)
  Cf_body <- skin_friction_cf(vs, aero$rocket_length, M)
  Cf_fin  <- skin_friction_cf(vs, aero$cbar,          M)
  Cd_fric <- (Cf_body * (aero$Awet_nose + (1 + 2/aero$fB) * aero$Awet_body) +
                Cf_fin  *  aero$Awet_fins) / aero$Aref     # eq. 3.85
  Cd_base <- live_base_drag(M)                            # eq. 3.94
  Cd_lug  <- aero$lug_k * cd_stag(M) * aero$A_lug / aero$Aref  # eq. 3.95
  
  # cd_scale is 1 in normal flight; Monte Carlo perturbs it. It multiplies the
  # WHOLE zero-AoA drag coefficient (including base drag), not just the
  # wetted-area terms, so the sampled drag uncertainty is the one the user asked for.
  list(CNa = CNa, CP = CP, Xcp_fin = Xcp_fin,
       CD0 = (Cd_pres + Cd_fric + Cd_base + Cd_lug) * aero$cd_scale)
}

# =============================================================================
# FLIGHT SIMULATION -- RK4, pitch-plane rigid body
# =============================================================================
# State vector y (13):
#   1:3   position   (x, y, z)      world, z up
#   4:6   velocity   (vx, vy, vz)   world
#   7:9   body axis unit vector u   world  (points out the nose)
#  10:12  angular velocity omega    world  (rad/s, perpendicular to u)
#  13     propellant mass remaining (kg)
#
# Sign convention for the normal force, which is what makes weathercocking
# come out the right way round:
#   w_hat = unit transverse component of the airflow direction, perpendicular
#           to the body axis.  The normal force is N * (-w_hat) -- it pushes
#           toward the side the nose is pointing (air strikes the windward
#           flank).  Applied at the CP, which lies (CP - CG) AFT of the CG,
#           this gives torque = N*(CP - CG) * (u x w_hat), which rotates u
#           TOWARD the relative wind.  A crosswind therefore turns the nose
#           INTO the wind (upwind weathercock), and a rocket whose CP is
#           forward of the CG diverges instead of correcting.
# =============================================================================

rocket_deriv <- function(t, y, ctx) {
  pos <- y[1:3]; vel <- y[4:6]; u <- y[7:9]; om <- y[10:12]
  m_prop <- max(y[13], 0)
  
  u <- unitv(u)
  z <- pos[3]
  
  m  <- ctx$eff_dry_mass + m_prop
  Ft <- ctx$thrust(t)
  
  # live CG and pitch inertia (uniform rod + parallel axis, Niskanen sec. 4.2.3)
  cg   <- (ctx$eff_dry_mass * ctx$cg_dry + m_prop * ctx$cg_motor) / m
  L    <- ctx$aero$rocket_length
  Ilong <- max(m * L^2 / 12 + m * (cg - L/2)^2, 1e-8)
  
  atm   <- isa_atm(z)
  rho   <- atm$rho
  a_snd <- atm$a
  
  vrel <- vel - ctx$wind
  vrm_raw <- vnorm(vrel)
  vrm  <- max(vrm_raw, 1e-6)
  vhat <- vrel / vrm
  M    <- vrm / a_snd
  q    <- 0.5 * rho * vrm^2
  
  dom  <- c(0, 0, 0)
  du   <- c(0, 0, 0)
  
  if (vrm_raw < 1e-3) {
    # No meaningful airflow: the airflow DIRECTION is undefined and the dynamic
    # pressure is negligible anyway (q < 1e-6 Pa). Coast on gravity + thrust.
    Faero <- c(0, 0, 0)
    
  } else if (ctx$chute_open) {
    # 3-DOF descent under canopy: all drag from the recovery device
    Faero <- -Cd_chute * ctx$chute_area * q * vhat
    Ft    <- 0
    
  } else if (ctx$on_rail) {
    # The rail holds the rocket rigidly along its axis, so there is no angle of
    # attack and no corrective moment yet. Only the AXIAL component of the
    # airflow produces force. (Without this, a stationary rocket in a crosswind
    # would register a spurious 90 deg AoA at t = 0.)
    v_ax  <- sum(vrel * u)
    M_ax  <- abs(v_ax) / a_snd
    co    <- aero_coeffs(ctx$aero, abs(v_ax), M_ax)
    Faero <- -co$CD0 * 0.5 * rho * v_ax^2 * ctx$aero$Aref * sign(v_ax) * u
    
  } else {
    co  <- aero_coeffs(ctx$aero, vrm, M)
    
    cosA  <- max(min(sum(u * vhat), 1), -1)
    alpha <- acos(cosA)                                  # 0 .. pi
    
    # axial force: CD0 scaled for angle of attack (sec. 3.4.7), along the body
    CA    <- co$CD0 * aoa_drag_factor(alpha)
    F_ax  <- -CA * q * ctx$aero$Aref * sign(cosA) * u
    
    # normal force: perpendicular to the body axis
    trans <- vhat - cosA * u
    tn    <- vnorm(trans)
    if (tn > 1e-9) {
      what <- trans / tn
      CN   <- co$CNa * sin(alpha)      # = CNa*alpha for small alpha; saturates
      F_n  <- -CN * q * ctx$aero$Aref * what
      torque <- CN * q * ctx$aero$Aref * (co$CP - cg) * cross3(u, what)
    } else {
      F_n <- c(0, 0, 0); torque <- c(0, 0, 0)
    }
    Faero <- F_ax + F_n
    
    # ---- pitch damping (eqs. 3.58-3.60) ----
    omn <- vnorm(om)
    if (omn > 1e-6) {
      l_f <- max(cg, 0); l_a <- max(L - cg, 0)
      Md_body <- 0.275 * rho * ctx$aero$bt_radius * (l_f^4 + l_a^4) * omn^2
      xi      <- abs(ctx$aero$Xcp_fin0 - cg)
      Md_fin  <- 0.3 * rho * min(ctx$aero$fin_count, 4) *
        ctx$aero$Afin_one * xi^3 * omn^2
      Md <- Md_body + Md_fin
      # damping can never reverse the rotation within one step
      Md <- min(Md, 0.5 * omn * Ilong / max(ctx$dt, 1e-6))
      torque <- torque - Md * (om / omn)
    }
    
    dom <- torque / Ilong
    du  <- cross3(om, u)
  }
  
  Fgrav <- c(0, 0, -m * g0)
  acc   <- (Faero + Fgrav + Ft * u) / m
  
  if (ctx$on_rail) {
    # The rail takes the lateral reaction: only motion along the rail survives,
    # and the rocket cannot rotate or slide backwards down the rail.
    a_par <- sum(acc * ctx$lvec)
    if (a_par < 0 && sum(vel * ctx$lvec) <= 0) a_par <- 0
    acc <- a_par * ctx$lvec
    dom <- c(0, 0, 0)
    du  <- c(0, 0, 0)
  }
  
  dmp <- if (t <= ctx$burn_time && m_prop > 0) -Ft / ctx$v_exhaust else 0
  
  c(vel, acc, du, dom, dmp)
}

flight_simulation_3d <- function(thrust_curve, prop_mass, dry_mass,
                                 casing_mass,
                                 bt_diameter, chute_diameter, chute_delay,
                                 aero, CNa_total, CP, precision,
                                 wind_speed_ref, wind_dir_deg,
                                 cg_dry_m, nose_length, body_length,
                                 motor_length_m, rail_length_m,
                                 launch_bearing_deg, launch_angle_deg,
                                 cd_scale    = 1.0,
                                 landing_only = FALSE,
                                 wind_turbulence_intensity = 15,
                                 gust_duration = 2) {
  
  if (is.null(thrust_curve) || nrow(thrust_curve) == 0) return(NULL)
  
  thrust    <- approxfun(thrust_curve$time, thrust_curve$thrust,
                         yleft = 0, yright = 0)
  burn_time <- max(thrust_curve$time)
  
  eff_dry_mass <- dry_mass + casing_mass
  
  total_imp <- integrate(thrust, min(thrust_curve$time),
                         max(thrust_curve$time))$value
  v_exhaust <- total_imp / prop_mass
  if (!is.finite(v_exhaust) || v_exhaust <= 0) return(NULL)
  
  chute_area <- pi * (chute_diameter / 2)^2
  cg_motor   <- nose_length + body_length - motor_length_m / 2
  
  # Monte Carlo drag perturbation: scales the entire CD0 inside aero_coeffs()
  aero_s <- aero
  aero_s$cd_scale <- cd_scale
  
  # launch rail direction
  ar <- launch_angle_deg   * pi / 180
  br <- launch_bearing_deg * pi / 180
  lvec <- c(sin(ar) * sin(br), sin(ar) * cos(br), cos(ar))
  lvec <- unitv(lvec)
  
  # ---- initial state: on the rail, aligned with it, no rotation ----
  y <- c(0, 0, 0,          # position
         0, 0, 0,          # velocity
         lvec,             # body axis
         0, 0, 0,          # angular velocity
         prop_mass)        # propellant
  
  t          <- 0
  dt         <- precision
  t_apogee   <- NA_real_
  t_eject    <- burn_time + chute_delay     # ejection: burnout + motor delay
  chute_open <- FALSE
  on_rail    <- TRUE
  rail_exit_v <- NA_real_
  deploy_v    <- NA_real_
  deploy_alt  <- NA_real_
  apogee_z    <- 0
  
  # ---- wind: Ornstein-Uhlenbeck turbulence about a sheared mean ----
  wd_rad  <- wind_dir_deg * pi / 180
  alpha_w <- 1 / max(gust_duration, 1e-3)
  wind    <- c(-wind_speed_ref * sin(wd_rad), -wind_speed_ref * cos(wd_rad), 0)
  
  max_steps <- ceiling(1200 / max(precision, 0.001))
  if (!landing_only) {
    out <- data.frame(
      time = numeric(max_steps), x = numeric(max_steps), y = numeric(max_steps),
      altitude = numeric(max_steps), velocity = numeric(max_steps),
      vx = numeric(max_steps), vy = numeric(max_steps), vz = numeric(max_steps),
      mach = numeric(max_steps), aoa = numeric(max_steps),
      stability_margin = numeric(max_steps), phase = integer(max_steps)
    )
    i <- 1L
  }
  
  repeat {
    z <- y[3]
    
    # --- wind update, once per step (held constant across the RK4 stages) ---
    z_ref  <- 10.0
    shear  <- (max(z, z_ref) / z_ref)^0.14          # 1/7 power law
    mu     <- c(-wind_speed_ref * shear * sin(wd_rad),
                -wind_speed_ref * shear * cos(wd_rad), 0)
    # Stationary sd of an OU process is sigma/sqrt(2*alpha). Scaling the
    # diffusion by sqrt(2*alpha) makes the realized turbulence intensity equal
    # the slider value for ANY gust duration (Niskanen eq. 4.7: I_u = sd/U).
    sigma_u <- 0.01 * wind_turbulence_intensity * wind_speed_ref * shear
    diff_c  <- sigma_u * sqrt(2 * alpha_w)
    wind[1] <- wind[1] + alpha_w*(mu[1] - wind[1])*dt + diff_c*sqrt(dt)*rnorm(1)
    wind[2] <- wind[2] + alpha_w*(mu[2] - wind[2])*dt + diff_c*sqrt(dt)*rnorm(1)
    wind[3] <- 0
    
    # --- events ---
    s_rail <- sum(y[1:3] * lvec)
    if (on_rail && s_rail >= rail_length_m) {
      on_rail <- FALSE
      rail_exit_v <- vnorm(y[4:6])
    }
    if (!chute_open && t >= t_eject) {
      chute_open <- TRUE
      deploy_v   <- vnorm(y[4:6])
      deploy_alt <- z
    }
    
    ctx <- list(aero = aero_s, thrust = thrust, burn_time = burn_time,
                v_exhaust = v_exhaust, eff_dry_mass = eff_dry_mass,
                cg_dry = cg_dry_m, cg_motor = cg_motor,
                chute_area = chute_area, chute_open = chute_open,
                on_rail = on_rail, lvec = lvec, wind = wind, dt = dt)
    
    # --- record ---
    if (!landing_only) {
      vrel  <- y[4:6] - wind
      vrm   <- max(vnorm(vrel), 1e-6)
      atm   <- isa_atm(z)
      uu    <- unitv(y[7:9])
      cosA  <- max(min(sum(uu * (vrel/vrm)), 1), -1)
      m_now <- eff_dry_mass + max(y[13], 0)
      cg_n  <- (eff_dry_mass * cg_dry_m + max(y[13],0) * cg_motor) / m_now
      Mn    <- vrm / atm$a
      cp_n  <- aero_coeffs(aero_s, vrm, Mn)$CP
      
      out$time[i]     <- t
      out$x[i]        <- y[1]; out$y[i] <- y[2]; out$altitude[i] <- z
      out$velocity[i] <- vnorm(y[4:6])
      out$vx[i]       <- y[4]; out$vy[i] <- y[5]; out$vz[i] <- y[6]
      out$mach[i]     <- Mn
      # AoA is undefined while the rail constrains the rocket (at t = 0 the only
      # relative airflow is the wind, which would read as 90 deg), and it is
      # equally undefined when the airspeed is ~0 -- which happens at apogee in
      # dead-still air. Report 0 in both cases; the aero forces there are nil.
      out$aoa[i]      <- if (on_rail || vrm < 1.0) 0 else acos(cosA) * 180 / pi
      out$stability_margin[i] <- (cp_n - cg_n) / bt_diameter
      out$phase[i] <- if (t <= burn_time) 1L else if (!chute_open) 2L else 3L
      i <- i + 1L
    }
    
    # --- RK4 (eqs. 4.20-4.21) ---
    k1 <- rocket_deriv(t,        y,                 ctx)
    k2 <- rocket_deriv(t + dt/2, y + k1 * (dt/2),   ctx)
    k3 <- rocket_deriv(t + dt/2, y + k2 * (dt/2),   ctx)
    k4 <- rocket_deriv(t + dt,   y + k3 * dt,       ctx)
    y  <- y + (dt/6) * (k1 + 2*k2 + 2*k3 + k4)
    t  <- t + dt
    
    # --- housekeeping: renormalize the body axis, keep omega perpendicular ---
    y[7:9] <- unitv(y[7:9])
    y[10:12] <- y[10:12] - sum(y[10:12] * y[7:9]) * y[7:9]
    y[13] <- max(y[13], 0)
    if (on_rail && sum(y[4:6] * lvec) < 0) y[4:6] <- c(0, 0, 0)
    
    if (y[3] > apogee_z) apogee_z <- y[3]
    if (is.na(t_apogee) && y[3] > 1 && y[6] < 0) t_apogee <- t
    
    # --- termination ---
    if (!is.na(t_apogee) && y[3] <= 0) {
      if (landing_only) {
        # linear interpolation back to the ground plane
        frac <- if (abs(y[6] * dt) > 1e-9) y[3] / (y[6] * dt) else 0
        return(data.frame(x = y[1] - y[4] * dt * frac,
                          y = y[2] - y[5] * dt * frac))
      }
      break
    }
    if (t > 1200 || (!landing_only && i >= max_steps)) break
    if (!landing_only && any(!is.finite(y))) break
    if (landing_only  && any(!is.finite(y))) return(data.frame(x = 0, y = 0))
  }
  
  if (landing_only) return(data.frame(x = y[1], y = y[2]))
  
  result <- out[1:(i - 1), ]
  attr(result, "rail_exit_ms") <- rail_exit_v
  attr(result, "t_eject")      <- t_eject
  attr(result, "t_apogee")     <- t_apogee
  attr(result, "deploy_v")     <- deploy_v
  attr(result, "deploy_alt")   <- deploy_alt
  attr(result, "burn_time")    <- burn_time
  result
}

# Rows from liftoff to apogee. Angle of attack is only meaningful as a stability
# diagnostic during the ascent: EVERY rocket, stable or not, swings through a
# large angle of attack as it arcs over at apogee, so measuring peak AoA over
# the whole coast phase would flag every single flight.
ascent_rows <- function(r) {
  ta <- attr(r, "t_apogee")
  if (is.null(ta) || is.na(ta)) r$phase == 1L else r$time <= ta
}
ascent_max_aoa <- function(r) {
  idx <- ascent_rows(r)
  if (!any(idx)) return(0)
  max(r$aoa[idx])
}

# =============================================================================
# MOTOR FILE PARSING (RASP .eng)
# =============================================================================
parse_thrust_input <- function(motor_file, engine_choice) {
  read_eng <- function(lines) {
    lines <- lines[!grepl("^;", lines)]
    lines <- lines[nzchar(trimws(lines))]
    hdr   <- strsplit(trimws(lines[1]), "\\s+")[[1]]
    # RASP header: name diam(mm) length(mm) delays prop_mass(kg) total_mass(kg) mfr
    pm       <- as.numeric(hdr[5])
    total_m  <- as.numeric(hdr[6])
    casing_m <- max(total_m - pm, 0)
    
    # field 4 is the delay list, e.g. "3-5-7", or "P"/"0" for plugged
    delays <- suppressWarnings(as.numeric(strsplit(hdr[4], "-")[[1]]))
    delays <- delays[!is.na(delays) & delays > 0]
    
    pairs <- lapply(lines[-1], function(l) as.numeric(strsplit(trimws(l), "\\s+")[[1]]))
    pairs <- Filter(function(p) length(p) >= 2 && !anyNA(p), pairs)
    list(
      thrust_curve   = data.frame(time = sapply(pairs, `[`, 1),
                                  thrust = sapply(pairs, `[`, 2)),
      prop_mass      = pm,
      casing_mass    = casing_m,
      delays         = delays,
      motor_length_m = as.numeric(hdr[3]) / 1000,
      motor_diam_m   = as.numeric(hdr[2]) / 1000
    )
  }
  if (!is.null(motor_file)) {
    tryCatch(read_eng(readLines(motor_file$datapath, warn=FALSE)), error=function(e) NULL)
  } else {
    if (is.null(engine_choice) || engine_choice == "") return(NULL)
    path <- file.path("data", paste0(engine_choice, ".eng"))
    if (!file.exists(path)) return(NULL)
    tryCatch(read_eng(readLines(path, warn=FALSE)), error=function(e) NULL)
  }
}

# =============================================================================
# OPENROCKET .ork IMPORT (best effort -- geometry only)
# .ork is a zip containing rocket.ork, an XML document. All lengths in metres.
# =============================================================================
parse_ork <- function(path) {
  tmp <- file.path(tempdir(), paste0("ork_", as.integer(runif(1, 1, 1e8))))
  dir.create(tmp, showWarnings = FALSE)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  
  files <- tryCatch(utils::unzip(path, exdir = tmp), error = function(e) character(0))
  xmlf  <- files[grepl("\\.ork$|\\.xml$", files)]
  doc   <- if (length(xmlf) > 0) {
    tryCatch(xml2::read_xml(xmlf[1]), error = function(e) NULL)
  } else {
    tryCatch(xml2::read_xml(path), error = function(e) NULL)   # uncompressed .ork
  }
  if (is.null(doc)) stop("Could not read .ork (not a valid OpenRocket file)")
  
  num <- function(node, tag) {
    if (is.na(node) || length(node) == 0) return(NA_real_)
    v <- xml2::xml_find_first(node, paste0("./", tag))
    if (length(v) == 0 || is.na(v)) return(NA_real_)
    suppressWarnings(as.numeric(xml2::xml_text(v)))
  }
  txt <- function(node, tag) {
    if (is.na(node) || length(node) == 0) return(NA_character_)
    v <- xml2::xml_find_first(node, paste0("./", tag))
    if (length(v) == 0 || is.na(v)) return(NA_character_)
    xml2::xml_text(v)
  }
  
  nose <- xml2::xml_find_first(doc, "//nosecone")
  tube <- xml2::xml_find_first(doc, "//bodytube")
  fins <- xml2::xml_find_first(doc, "//trapezoidfinset")
  
  shape_raw <- tolower(txt(nose, "shape"))
  nose_type <- if (is.na(shape_raw)) "ogive"
  else if (grepl("cone", shape_raw))      "conical"
  else if (grepl("parab|ellip", shape_raw)) "parabolic"
  else                                     "ogive"
  
  nose_len <- num(nose, "length")
  aftrad   <- num(nose, "aftradius")
  body_len <- num(tube, "length")
  tube_rad <- num(tube, "radius")
  radius   <- if (!is.na(tube_rad)) tube_rad else aftrad
  
  res <- list(
    nose_type   = nose_type,
    nose_length = nose_len,
    body_length = body_len,
    diameter    = if (!is.na(radius)) 2 * radius else NA_real_,
    fin_count   = num(fins, "fincount"),
    fin_root    = num(fins, "rootchord"),
    fin_tip     = num(fins, "tipchord"),
    fin_span    = num(fins, "height"),
    fin_sweep   = num(fins, "sweeplength")
  )
  # Fin axial position: OpenRocket stores it relative to the tube; if the
  # fin set is bottom-referenced, the root LE sits root-chord forward of the aft end.
  if (!is.na(res$nose_length) && !is.na(res$body_length) && !is.na(res$fin_root)) {
    res$fin_pos <- res$nose_length + res$body_length - res$fin_root
  }
  # Mass / CG overrides if the designer set them
  om <- xml2::xml_find_first(doc, "//overridemass")
  oc <- xml2::xml_find_first(doc, "//overridecg")
  if (length(om) > 0 && !is.na(om))
    res$dry_mass <- suppressWarnings(as.numeric(xml2::xml_text(om)))
  if (length(oc) > 0 && !is.na(oc))
    res$cg_measured <- suppressWarnings(as.numeric(xml2::xml_text(oc)))
  
  res[!vapply(res, function(v) length(v) == 0 || all(is.na(v)), logical(1))]
}

# ---- OpenRocket CSV export (validation overlay) ------------------------------
parse_or_csv <- function(path) {
  df <- tryCatch(
    utils::read.csv(path, comment.char = "#", check.names = FALSE,
                    stringsAsFactors = FALSE),
    error = function(e) NULL)
  if (is.null(df) || ncol(df) < 2) return(NULL)
  nm  <- tolower(names(df))
  tcol <- which(grepl("time", nm))[1]
  acol <- which(grepl("altitude|apogee|height", nm))[1]
  if (is.na(tcol) || is.na(acol)) return(NULL)
  out <- data.frame(time = suppressWarnings(as.numeric(df[[tcol]])),
                    altitude = suppressWarnings(as.numeric(df[[acol]])))
  out <- out[is.finite(out$time) & is.finite(out$altitude), ]
  if (nrow(out) == 0) return(NULL)
  out
}

# ---- CSS (unchanged) --------------------------------------------------------
css <- "
@import url('https://fonts.googleapis.com/css2?family=Lexend:wght@300;400;500;600;700;800&display=swap');

:root {
  --c1:    #ff4e50;
  --c2:    #fc913a;
  --c3:    #f9d62e;
  --c4:    #eae374;
  --c5:    #e2f4c7;
  --dark:  #1a1008;
  --panel: rgba(26,16,8,0.82);
  --border:#fc913a44;
  --text:  #fff8f0;
  --dim:   #f9d62ecc;
  --sans:  'Lexend', sans-serif;
}

*,*::before,*::after { font-family:'Lexend',sans-serif!important; }

html {
  min-height:100%;
  background:
    linear-gradient(180deg,
      #ff4e50 0%,   #ff4e50 16%,
      #fc913a 16%,  #fc913a 32%,
      #f9d62e 32%,  #f9d62e 48%,
      #eae374 48%,  #eae374 64%,
      #e2f4c7 64%,  #e2f4c7 80%,
      #b8dba0 80%,  #b8dba0 100%) fixed!important;
}
body {
  background:transparent!important;
  color:var(--text)!important;
  min-height:100vh;
}
body::before {
  content:'';position:fixed;inset:0;
  background:repeating-linear-gradient(180deg,transparent 0px,transparent calc(16vh - 4px),
    rgba(0,0,0,0.22) calc(16vh - 4px),rgba(255,255,255,0.06) calc(16vh),transparent calc(16vh + 1px));
  pointer-events:none;z-index:0;
}
body::after {
  content:'';position:fixed;inset:0;background:rgba(10,6,2,0.55);pointer-events:none;z-index:1;
}
.tab-content,.navbar,.shiny-plot-output,.card,.well,.leaflet-container,
pre,.shiny-verbatim-output,.stab-box,.run-table,.irs--shiny,
.form-control,input,select,.nav-pills,.sidebarPanel,.mainPanel,.container-fluid {
  position:relative;z-index:2;
}

/* NAVBAR */
.navbar {
  background:linear-gradient(180deg,rgba(26,16,8,0.97) 0%,rgba(20,10,4,0.99) 100%)!important;
  border-bottom:2px solid var(--c2)!important;
  box-shadow:0 4px 24px rgba(255,78,80,0.25),0 2px 8px rgba(0,0,0,0.7)!important;
  padding:4px 32px!important;z-index:100!important;position:relative!important;min-height:58px!important;
}
.navbar-brand { display:flex!important;align-items:center!important;padding:8px 0!important;gap:0!important; }
.nav-link {
  font-size:0.73rem!important;color:#f9d62eaa!important;padding:16px 20px!important;
  text-transform:uppercase!important;letter-spacing:1.4px!important;font-weight:600!important;
  border-bottom:3px solid transparent!important;transition:all .2s!important;
}
.nav-link:hover { color:var(--c3)!important; }
.nav-link.active {
  color:var(--c1)!important;border-bottom-color:var(--c1)!important;
  background:transparent!important;text-shadow:0 0 12px rgba(255,78,80,0.5)!important;
}

/* BUTTONS */
.btn {
  font-size:0.73rem!important;text-transform:uppercase!important;letter-spacing:1.1px!important;
  font-weight:700!important;border-radius:7px!important;padding:10px 26px!important;
  transition:all .08s ease!important;cursor:pointer!important;
}
.btn-primary {
  background:linear-gradient(180deg,#ff6b5b 0%,#ff4e50 50%,#d93535 100%)!important;
  border:1px solid #b02020!important;color:#fff!important;
  box-shadow:0 1px 0 rgba(255,255,255,0.3) inset,0 -2px 0 rgba(0,0,0,0.3) inset,
    0 4px 10px rgba(255,78,80,0.5),0 2px 4px rgba(0,0,0,0.5)!important;
  text-shadow:0 -1px 0 rgba(0,0,0,0.4)!important;
}
.btn-primary:hover {
  background:linear-gradient(180deg,#ff8070 0%,#ff6b5b 50%,#ff4e50 100%)!important;
  box-shadow:0 1px 0 rgba(255,255,255,0.35) inset,0 -2px 0 rgba(0,0,0,0.25) inset,
    0 6px 16px rgba(255,78,80,0.6),0 2px 6px rgba(0,0,0,0.4)!important;
}
.btn-primary:active {
  background:linear-gradient(180deg,#d93535 0%,#b02020 100%)!important;
  box-shadow:0 2px 4px rgba(0,0,0,0.6) inset!important;transform:translateY(1px)!important;
}
.btn-warning {
  background:linear-gradient(180deg,#ffaa55 0%,#fc913a 50%,#d97020 100%)!important;
  border:1px solid #b05810!important;color:#fff!important;
  box-shadow:0 1px 0 rgba(255,255,255,0.3) inset,0 -2px 0 rgba(0,0,0,0.3) inset,
    0 4px 10px rgba(252,145,58,0.5),0 2px 4px rgba(0,0,0,0.5)!important;
  text-shadow:0 -1px 0 rgba(0,0,0,0.35)!important;
}
.btn-warning:hover {
  background:linear-gradient(180deg,#ffc070 0%,#ffaa55 50%,#fc913a 100%)!important;
  box-shadow:0 1px 0 rgba(255,255,255,0.35) inset,0 -2px 0 rgba(0,0,0,0.25) inset,
    0 6px 16px rgba(252,145,58,0.55),0 2px 6px rgba(0,0,0,0.4)!important;
}
.btn-warning:active {
  background:linear-gradient(180deg,#d97020 0%,#b05810 100%)!important;
  box-shadow:0 2px 4px rgba(0,0,0,0.6) inset!important;transform:translateY(1px)!important;
}
.btn-default,.btn-secondary {
  background:linear-gradient(180deg,#3a2510 0%,#251508 100%)!important;
  border:1px solid #fc913a66!important;color:#f9d62e!important;
  box-shadow:0 1px 0 rgba(255,255,255,0.08) inset,0 -1px 0 rgba(0,0,0,0.4) inset,
    0 2px 6px rgba(0,0,0,0.5)!important;
}
.btn-default:hover,.btn-secondary:hover {
  background:linear-gradient(180deg,#4a3018 0%,#3a2510 100%)!important;color:var(--c3)!important;
}
.btn-file {
  background:linear-gradient(180deg,#3a2510,#251508)!important;
  border:1px solid #fc913a55!important;color:#f9d62e!important;box-shadow:0 2px 6px rgba(0,0,0,0.5)!important;
}

/* INPUTS */
label,.form-label {
  font-size:0.68rem!important;color:var(--c4)!important;text-transform:uppercase!important;
  letter-spacing:.08em!important;font-weight:600!important;
}
.form-control,.form-select,input[type=number],select {
  background:rgba(20,10,4,0.85)!important;border:1px solid #fc913a55!important;
  color:var(--text)!important;border-radius:5px!important;padding:6px 10px!important;
  font-size:0.84rem!important;
  box-shadow:0 2px 5px rgba(0,0,0,0.5) inset,0 1px 0 rgba(255,255,255,0.04)!important;
  transition:border-color .15s,box-shadow .15s!important;
}
.form-control:focus,input:focus,select:focus {
  border-color:var(--c2)!important;
  box-shadow:0 2px 5px rgba(0,0,0,0.5) inset,0 0 0 3px rgba(252,145,58,0.25)!important;
  outline:none!important;
}
select option { background:#1a0f05;color:var(--text); }

/* SLIDERS */
.irs--shiny .irs-bar { background:linear-gradient(90deg,var(--c1),var(--c2))!important;border-color:var(--c1)!important; }
.irs--shiny .irs-handle { border-color:var(--c2)!important;background:#fff8f0!important;box-shadow:0 2px 6px rgba(0,0,0,0.5)!important; }
.irs--shiny .irs-single { background:var(--c2)!important;color:#1a1008!important;font-size:0.65rem!important;font-weight:700!important; }
.irs--shiny .irs-line { background:#2a1508!important;box-shadow:0 1px 4px rgba(0,0,0,0.5) inset!important; }
.irs--shiny .irs-grid-text,.irs--shiny .irs-min,.irs--shiny .irs-max { color:var(--c4)!important;background:transparent!important; }
.irs--shiny .irs-grid-pol { background:var(--c4)!important; }

/* CARDS */
.card {
  background:linear-gradient(160deg,rgba(30,15,5,0.92) 0%,rgba(18,9,3,0.95) 100%)!important;
  border:1px solid #fc913a44!important;border-radius:8px!important;
  box-shadow:0 6px 20px rgba(0,0,0,0.6),0 1px 0 rgba(255,200,100,0.06) inset!important;
}
.card-body { padding:16px!important; }
.well {
  background:rgba(18,9,3,0.88)!important;border:1px solid #fc913a44!important;
  border-radius:8px!important;box-shadow:0 3px 10px rgba(0,0,0,0.5) inset!important;
}

/* NAV PILLS */
.nav-pills .nav-link {
  font-size:0.67rem!important;color:var(--c4)!important;border-radius:5px!important;
  padding:5px 13px!important;text-transform:uppercase!important;letter-spacing:.08em!important;
  background:linear-gradient(180deg,#3a2010 0%,#1e0e04 100%)!important;
  border:1px solid #fc913a44!important;margin:2px!important;font-weight:600!important;
  box-shadow:0 2px 5px rgba(0,0,0,0.4),0 1px 0 rgba(255,200,80,0.07) inset!important;
  transition:all .12s!important;
}
.nav-pills .nav-link:hover { background:linear-gradient(180deg,#4a2a14 0%,#2e1608 100%)!important;color:var(--c3)!important; }
.nav-pills .nav-link.active {
  background:linear-gradient(180deg,#ff6b5b 0%,#ff4e50 50%,#d93535 100%)!important;
  border-color:#b02020!important;color:#fff!important;
  box-shadow:0 3px 10px rgba(255,78,80,0.45),0 1px 0 rgba(255,255,255,0.2) inset!important;
}

p { font-size:0.7rem;color:var(--c4);margin:2px 0 8px; }
.tab-content { padding:20px 24px; }

pre,.shiny-verbatim-output {
  background:rgba(8,4,0,0.92)!important;color:#f9d62e!important;
  border:1px solid #fc913a44!important;border-radius:6px!important;
  font-family:'Courier New',monospace!important;font-size:0.78rem!important;
  padding:14px!important;line-height:1.9!important;box-shadow:0 3px 10px rgba(0,0,0,0.6) inset!important;
}
.shiny-plot-output { border:1px solid #fc913a44;border-radius:6px;box-shadow:0 6px 16px rgba(0,0,0,0.5); }
.leaflet-container { border-radius:6px;border:1px solid #fc913a55;box-shadow:0 6px 16px rgba(0,0,0,0.5); }

/* HOME */
.home-wrap { display:flex;flex-direction:column;justify-content:center;align-items:center;min-height:80vh;text-align:center;position:relative; }
.home-graphic {
  position:absolute;width:520px;height:520px;border-radius:50%;
  background:radial-gradient(ellipse at 50% 60%,rgba(249,214,46,0.18) 0%,rgba(252,145,58,0.14) 35%,rgba(255,78,80,0.10) 60%,transparent 75%);
  box-shadow:0 0 80px 20px rgba(249,214,46,0.12),0 0 160px 60px rgba(252,145,58,0.08),0 0 260px 100px rgba(255,78,80,0.05);
  top:50%;left:50%;transform:translate(-50%,-50%);pointer-events:none;z-index:0;
}
.home-horizon {
  position:absolute;bottom:28%;left:0;right:0;height:1px;
  background:linear-gradient(90deg,transparent 0%,rgba(249,214,46,0.35) 20%,rgba(252,145,58,0.6) 50%,rgba(249,214,46,0.35) 80%,transparent 100%);
  box-shadow:0 0 12px 2px rgba(252,145,58,0.3);pointer-events:none;
}
.home-mountains {
  position:absolute;bottom:28%;left:0;right:0;height:80px;
  clip-path:polygon(0% 100%,6% 55%,12% 75%,20% 30%,28% 60%,36% 20%,44% 50%,52% 35%,60% 55%,68% 10%,76% 45%,84% 30%,92% 50%,100% 40%,100% 100%);
  background:#1a0a02;opacity:0.85;pointer-events:none;z-index:1;
}
.home-content { position:relative;z-index:2; }
.title-wrap { display:flex;align-items:center;justify-content:center; }
.home-title {
  font-size:5rem;font-weight:800;letter-spacing:3px;margin:0;
  background:linear-gradient(135deg,#ff4e50 0%,#fc913a 35%,#f9d62e 65%,#e2f4c7 100%);
  -webkit-background-clip:text;-webkit-text-fill-color:transparent;
  text-rendering:geometricPrecision;filter:drop-shadow(0 0 20px rgba(252,145,58,0.4));
}
.cursor {
  display:inline-block;width:4px;height:5rem;
  background:linear-gradient(180deg,#f9d62e,#fc913a);margin-left:8px;vertical-align:middle;
  animation:blink .7s step-end infinite;box-shadow:0 0 10px rgba(249,214,46,0.6);border-radius:2px;
}
@keyframes blink { 50% { opacity:0; } }
.home-sub { font-size:0.85rem;letter-spacing:4px;margin-top:12px;text-transform:uppercase;font-weight:500;color:var(--c4);text-shadow:0 0 20px rgba(234,227,116,0.4); }

/* STABILITY BOX */
.stab-box {
  border-radius:6px;padding:12px 16px;margin-top:12px;font-size:0.78rem;line-height:2;
  font-weight:500;box-shadow:0 4px 14px rgba(0,0,0,0.5);
}
.stab-green  { color:#86efac;border:1px solid #22c55e;background:linear-gradient(160deg,rgba(13,40,24,0.9),rgba(19,34,24,0.92)); }
.stab-yellow { color:#fde68a;border:1px solid #f59e0b;background:linear-gradient(160deg,rgba(43,31,8,0.9),rgba(43,33,17,0.92)); }
.stab-red    { color:#fca5a5;border:1px solid #ef4444;background:linear-gradient(160deg,rgba(42,14,14,0.9),rgba(42,21,21,0.92)); }

/* RUN TABLE */
.run-table { width:100%;border-collapse:collapse;font-size:0.73rem;margin-top:10px;table-layout:fixed;word-break:break-word; }
.run-table th { color:var(--c4);text-transform:uppercase;letter-spacing:.07em;padding:6px 10px;border-bottom:2px solid #fc913a55;text-align:left;font-weight:700; }
.run-table td { color:var(--text);padding:6px 10px;border-bottom:1px solid rgba(252,145,58,0.15);font-weight:500;overflow:hidden;text-overflow:ellipsis;white-space:nowrap; }
.run-table tr:hover td { background:rgba(252,145,58,0.07); }

/* UNIT INPUT */
.unit-row-wrap { margin-bottom:10px; }
.unit-row { display:flex;gap:6px;align-items:center; }
.num-wrap { flex:2; }
.sel-wrap { flex:1; }
.unit-lbl { font-size:0.67rem;color:var(--c4);text-transform:uppercase;letter-spacing:.08em;font-weight:600;display:block;margin-bottom:3px; }

/* FOOTER */
.app-footer {
  position:relative;z-index:2;
  border-top:1px solid #fc913a33;
  padding:18px 32px;
  margin-top:40px;
  display:flex;
  justify-content:space-between;
  align-items:center;
  font-size:0.65rem;
  color:#eae374aa;
  letter-spacing:0.06em;
}
.app-footer a { color:#fc913a88;text-decoration:none; }
.app-footer a:hover { color:var(--c2); }
"

# ---- UI helper: paired numeric + unit selector ------------------------------
unit_input <- function(input_id, label, default_val, default_unit, choices) {
  div(class="unit-row-wrap",
      tags$label(class="unit-lbl", label),
      div(class="unit-row",
          div(class="num-wrap", numericInput(input_id, label=NULL, value=default_val, width="100%")),
          div(class="sel-wrap", selectInput(paste0(input_id,"_unit"), label=NULL,
                                            choices=choices, selected=default_unit, width="100%"))
      )
  )
}

# ---- UI ---------------------------------------------------------------------
ui <- tagList(
  tags$head(tags$style(HTML(css))),
  navbarPage(
    title = span(style=paste0(
      "font-size:1.4rem;font-weight:800;letter-spacing:3px;",
      "background:linear-gradient(90deg,#ff4e50,#fc913a,#f9d62e);",
      "-webkit-background-clip:text;-webkit-text-fill-color:transparent;",
      "vertical-align:middle;line-height:1;"
    ), "RRRocket 3D"),
    theme = bs_theme(version=5, bg="#f5f6f8", fg="#111928", primary="#1a56db"),
    
    tabPanel("Home",
             tags$script(HTML("
        $(document).on('shiny:sessioninitialized', function() {
          var text = 'RRRocket 3D';
          var el = document.getElementById('typed-title');
          var i = 0; el.textContent = '';
          function type() { if (i < text.length) { el.textContent += text[i++]; setTimeout(type, 100); } }
          setTimeout(type, 300);
        });
      ")),
             div(class="home-wrap",
                 div(class="title-wrap",
                     h1(class="home-title", span(id="typed-title")),
                     span(class="cursor")),
                 p(class="home-sub", "Model rocket flight simulator")
             )
    ),
    
    tabPanel("Setup",
             fluidRow(
               column(6,
                      navset_card_pill(
                        nav_panel("Rocket",
                                  unit_input("dry_mass_val", "Dry mass", 90,  "g",  unit_choices_mass),
                                  unit_input("diameter",     "Body tube diameter",  24,  "mm", unit_choices_length),
                                  unit_input("body_length",  "Body tube length",    300, "mm", unit_choices_length),
                                  unit_input("cg_measured",  "CG from nose tip (dry)", 220, "mm", unit_choices_length)
                        ),
                        nav_panel("Nosecone",
                                  selectInput("nose_type","Nosecone type", choices=c("ogive","conical","parabolic")),
                                  unit_input("nose_length","Nosecone length", 70, "mm", unit_choices_length)
                        ),
                        nav_panel("Fins",
                                  numericInput("fin_count","Number of fins", value=3, min=1, max=12),
                                  unit_input("fin_root",  "Root chord",      50, "mm", unit_choices_length),
                                  unit_input("fin_tip",   "Tip chord",       25, "mm", unit_choices_length),
                                  unit_input("fin_span",  "Semi-span", 30, "mm", unit_choices_length),
                                  unit_input("fin_sweep", "Sweep length", 20, "mm", unit_choices_length),
                                  unit_input("fin_pos",   "Root leading edge from nose tip", 320, "mm", unit_choices_length)
                        ),
                        nav_panel("Chute",
                                  unit_input("parachute_diameter","Chute diameter", 305, "mm", unit_choices_length)
                        ),
                        nav_panel("Wind",
                                  unit_input("wind_speed_val", "Wind speed at 10m", 3, "m/s", unit_choices_speed),
                                  p("Wind increases with altitude per the 1/7 power law."),
                                  sliderInput("wind_dir","Direction wind is coming from (deg CW from N)", min=0, max=360, value=270),
                                  sliderInput("wind_turbulence_intensity", "Wind turbulence intensity % (10-20 typical)",
                                              min = 5, max = 50, value = 15),
                                  sliderInput("gust_duration", "Mean gust duration (s)", min = 0.5, max = 5, value = 2, step = 0.5),
                                  p("Turbulence sd is held equal to intensity x wind speed regardless of gust duration.")
                        ),
                        nav_panel("Launch Site",
                                  unit_input("rail_length", "Rail length", 0.9, "m", unit_choices_length),
                                  numericInput("launch_angle","Launch angle from vertical (deg)", value=0, min=0, max=30),
                                  sliderInput("launch_bearing","Launch bearing (deg CW from N)", value=0, min=0, max=360),
                                  tags$hr(style="border-color:#fc913a33;"),
                                  unit_input("lug_length", "Launch lug length", 25,  "mm", unit_choices_length),
                                  unit_input("lug_od",     "Launch lug outer diameter",    4.0, "mm", unit_choices_length),
                                  unit_input("lug_id",     "Launch lug inner diameter",    3.2, "mm", unit_choices_length)
                        ),
                        nav_panel("Engine",
                                  fileInput("motor_file", NULL, accept=".eng", buttonLabel="Upload .eng"),
                                  selectInput("engine_choice","Or choose an engine:",
                                              choices=c("Select engine..."="","A8","A10","B4","B6","C6","C11",
                                                        "D12","E12","E16","G40"),
                                              selected="B6", size=6, selectize=FALSE),
                                  numericInput("parachute_delay","Ejection Delay (s)", value=4, min=0),
                                  uiOutput("delay_hint")
                        ),
                        nav_panel("Design",
                                  p("Save the current design, reload it later, or import geometry from an OpenRocket file."),
                                  downloadButton("save_design", "Save design (.json)", class="btn-default"),
                                  br(), br(),
                                  fileInput("load_design", "Load design (.json)", accept=".json"),
                                  fileInput("ork_file", "Import OpenRocket (.ork)", accept=".ork"),
                                  p("The .ork import reads geometry only (nose, body, trapezoidal fins). Masses and CG must still be entered or measured.")
                        )
                      ),
                      br(),
                      uiOutput("stability_indicator")
               ),
               column(6,
                      plotOutput("thrust_curve_plot", height="400px"),
                      br(),
                      plotOutput("fin_preview", height = "380px")
               )
             )
    ),
    
    tabPanel("Simulate",
             br(),
             fluidRow(
               column(3,
                      div(class="card", style="padding:16px;",
                          h6(style="color:var(--c3);text-transform:uppercase;letter-spacing:1px;font-size:0.72rem;margin-bottom:12px;",
                             "Settings"),
                          sliderInput("precision","Integration interval (s)", value=0.02, min=0.005, max=0.1, step=0.005),
                          p("RK4: 0.02 s recommended."),
                          div(style="margin:10px 0 6px;", tags$label(class="unit-lbl","Display units")),
                          radioButtons("units", label=NULL,
                                       choices=c("Metric"="metric","Imperial (ft)"="imperial"),
                                       selected="metric", inline=FALSE),
                          br(),
                          actionButton("run","> Simulate", class="btn-primary", style="width:100%;"),
                          br(), br(),
                          fileInput("or_csv", "Overlay OpenRocket CSV (SI)", accept=".csv")
                      ),
                      br(),
                      conditionalPanel("output.has_results",
                                       div(class="card", style="padding:16px;",
                                           h6(style="color:var(--c3);text-transform:uppercase;letter-spacing:1px;font-size:0.72rem;margin-bottom:10px;",
                                              "Flight summary"),
                                           verbatimTextOutput("summary")
                                       ),
                                       br(),
                                       uiOutput("safety_box")
                      )
               ),
               column(9,
                      conditionalPanel("output.has_results",
                                       fluidRow(
                                         column(6, plotOutput("altitude_plot", height="260px")),
                                         column(6, plotOutput("velocity_plot", height="260px"))
                                       ),
                                       br(),
                                       fluidRow(
                                         column(6, plotOutput("aoa_plot", height="220px")),
                                         column(6, plotOutput("stab_plot", height="220px"))
                                       ),
                                       br(),
                                       plotlyOutput("track_3d", height="420px")
                      ),
                      conditionalPanel("!output.has_results",
                                       div(style=paste0(
                                         "display:flex;align-items:center;justify-content:center;",
                                         "height:500px;color:#fc913a44;font-size:0.85rem;",
                                         "text-transform:uppercase;letter-spacing:2px;"
                                       ), ">  Press Simulate to run a flight")
                      )
               )
             ),
             br(),
             conditionalPanel("output.has_results",
                              div(class="card", style="padding:16px;overflow-x:auto;",
                                  h6(style="color:var(--c3);text-transform:uppercase;letter-spacing:1px;font-size:0.72rem;margin-bottom:10px;",
                                     "Run history"),
                                  uiOutput("run_history_table")
                              )
             )
    ),
    
    tabPanel("Monte Carlo",
             sidebarLayout(
               sidebarPanel(
                 numericInput("mc_runs","Monte Carlo runs", value=200, min=10, max=1000),
                 p("RK4 costs ~4 force evaluations per step; use 0.05 s here."),
                 sliderInput("chute_delay_std_dev",  "Ejection delay sd (s)",    value=0.5,min=0.1, max=5, step=0.1),
                 sliderInput("launch_angle_std_dev", "Launch angle sd (deg)",    value=5,  min=0.1, max=10),
                 sliderInput("wind_speed_std_dev",   "Wind speed sd (%)",        value=5,  min=1,   max=99),
                 sliderInput("wind_dir_std_dev",     "Wind direction sd (deg)",  value=30, min=1,   max=180),
                 sliderInput("dry_mass_std_dev",     "Dry mass sd (%)",          value=2,  min=0,   max=10),
                 sliderInput("cd_std_dev",           "Drag coefficient sd (%)",  value=10, min=0,   max=30),
                 sliderInput("prop_mass_std_dev",    "Propellant mass sd (%)",   value=2,  min=0,   max=10),
                 sliderInput("montecarlo_precision", "Integration interval (s)", value=0.05, min=0.01, max=0.2)
               ),
               mainPanel(
                 leafletOutput("map", height=600),
                 br(),
                 h6("Click map to draw landing site polygon. Click again to set launchpad position (pin)."),
                 br(),
                 actionButton("run_mc","> Run Monte Carlo", class="btn-warning"),
                 conditionalPanel("output.has_mc",
                                  div(class="card", style="padding:14px;margin-top:12px;",
                                      verbatimTextOutput("landing_pct"))
                 )
               )
             )
    ),
    
    tabPanel("About",
             div(style="max-width:700px;margin:40px auto;",
                 div(class="card", style="padding:30px;",
                     h3(style="color:var(--c2);font-weight:800;letter-spacing:1px;margin-bottom:4px;","RRRocket 3D"),
                     p(style="color:var(--c4);font-size:0.8rem;letter-spacing:2px;text-transform:uppercase;margin-bottom:24px;",
                       "Model rocket flight simulator"),
                     tags$hr(style="border-color:#fc913a33;margin-bottom:18px;"),
                     h6(style="color:var(--c3);text-transform:uppercase;letter-spacing:1px;font-size:0.72rem;","What it does"),
                     p(style="color:var(--text);font-size:0.85rem;line-height:1.8;",
                       "RRRocket 3D simulates model rocket flights and landings. My goal is to create a website that is both accurate and easy-to-use for hobbyists and professionals alike."),
                     br(),
                     h6(style="color:var(--c3);text-transform:uppercase;letter-spacing:1px;font-size:0.72rem;","Physics"),
                     p(style="color:var(--text);font-size:0.85rem;line-height:1.8;",
                       "Aerodynamics follow Barrowman (1967) and the OpenRocket technical documentation (Niskanen 2013). The rocket is split into nosecone, body tube and fins, and the drag of each is computed separately."),
                     br(),
                     p(style="color:var(--text);font-size:0.85rem;line-height:1.8;",
                       "Drag is decomposed into nose pressure drag (half-angle sine-squared), skin friction drag (turbulent flat-plate, eq. 3.78, switching to the roughness-limited value above the critical Reynolds number, eq. 3.80, with subsonic and supersonic compressibility corrections, eqs. 3.82-3.84), base drag (Hoerner, eq. 3.94) and launch-lug parasitic drag (eqs. 3.95-3.96, scaled by the stagnation pressure coefficient of eq. B.2). The zero-angle drag coefficient is scaled with angle of attack, rising to 1.3x at 17 degrees and falling to zero at 90 degrees (sec. 3.4.7)."),
                     br(),
                     p(style="color:var(--text);font-size:0.85rem;line-height:1.8;",
                       "The rocket is simulated as a rigid body in the pitch plane. The state carries the body axis and the angular velocity, so the normal force acting at the CP produces a real pitching moment about the live CG. Weathercocking, gravity turn, angle-of-attack drag and the divergence of an unstable rocket all emerge from this moment rather than from a tuned gain. Fin normal force uses eq. 3.40 with the Prandtl factor, and the fin CP marches aft above Mach 0.5 (eqs. 3.35-3.36). Pitch damping (eqs. 3.58-3.60) suppresses the wild oscillation that would otherwise follow apogee."),
                     br(),
                     p(style="color:var(--text);font-size:0.85rem;line-height:1.8;",
                       "The atmosphere follows the ISA troposphere model. Wind uses a 1/7 power law shear profile with an Ornstein-Uhlenbeck turbulence process whose stationary standard deviation equals the specified turbulence intensity times the local wind speed, independent of the gust duration. Flight is integrated with Runge-Kutta 4 (eqs. 4.20-4.21). Mass, CG, stability margin and all aerodynamic coefficients are recomputed at every stage."),
                     br(),
                     h6(style="color:var(--c3);text-transform:uppercase;letter-spacing:1px;font-size:0.72rem;","Recovery"),
                     p(style="color:var(--text);font-size:0.85rem;line-height:1.8;",
                       "The ejection charge fires at motor burnout plus the motor's delay, exactly as a real motor does, and NOT at apogee. The simulator reports how far from apogee the charge fires and how fast the rocket is moving when the chute opens, and warns when deployment would be fast enough to damage the airframe. Under canopy the drag coefficient is 0.80 over the canopy area (Hoerner)."),
                     br(),
                     h6(style="color:var(--c3);text-transform:uppercase;letter-spacing:1px;font-size:0.72rem;","Monte Carlo"),
                     p(style="color:var(--text);font-size:0.85rem;line-height:1.8;",
                       "Monte Carlo perturbs ejection delay, launch angle, wind speed and direction, dry mass, drag coefficient and propellant mass with Gaussian samples about the nominal values. Each run returns only the landing coordinate. Results are plotted on a satellite map and the proportion landing inside a user-drawn polygon is reported."),
                     br(),
                     h6(style="color:var(--c3);text-transform:uppercase;letter-spacing:1px;font-size:0.72rem;","Validation"),
                     p(style="color:var(--text);font-size:0.85rem;line-height:1.8;",
                       "Designs can be saved and reloaded as JSON, imported from OpenRocket .ork files, and an OpenRocket CSV export can be overlaid on the altitude plot to compare the two simulators directly."),
                     br(),
                     h6(style="color:var(--c3);text-transform:uppercase;letter-spacing:1px;font-size:0.72rem;","Engine data"),
                     p(style="color:var(--text);font-size:0.85rem;line-height:1.8;",
                       "Motors use the standard RASP .eng format; the delay list in the header populates the ejection delay. For anything larger than the built-in Estes curves, download the .eng file from ",
                       tags$a(href="https://www.thrustcurve.org", target="_blank",
                              style="color:var(--c2);", "thrustcurve.org"), "."),
                     br(),
                     tags$hr(style="border-color:#fc913a33;margin:8px 0 20px;"),
                     div(style="display:flex;gap:12px;",
                         tags$a(href="https://github.com/tatecommission/rrrocket", target="_blank",
                                class="btn btn-default","GitHub"),
                         tags$a(href="https://www.thrustcurve.org", target="_blank",
                                class="btn btn-default","ThrustCurve.org")
                     )
                 )
             )
    ),
    
    nav_spacer(),
    nav_item(
      tags$a(href="https://github.com/tatecommission/rrrocket", target="_blank",
             style=paste0(
               "font-size:0.65rem;font-weight:700;letter-spacing:1px;",
               "text-transform:uppercase;color:#eae374;text-decoration:none;",
               "border:1px solid #fc913a55;border-radius:5px;padding:4px 10px;",
               "background:rgba(252,145,58,0.1);",
               "box-shadow:0 1px 0 rgba(255,255,255,0.08) inset,0 2px 5px rgba(0,0,0,0.4);",
               "display:inline-block;"
             ), "<GitHub>"
      )
    ),
    tags$footer(class = "app-footer",
                span("\u00a9 2026 Tate Commission. All rights reserved."),
                span(
                  tags$a(href="https://github.com/tatecommission/rrrocket", target="_blank", "GitHub")
                )
    )
  )
)

# ---- plot theme -------------------------------------------------------------
theme_plot <- function() {
  theme_minimal(base_size=11) + theme(
    plot.background  = element_rect(fill="#120800", color=NA),
    panel.background = element_rect(fill="#120800", color=NA),
    panel.grid.major = element_line(color="#3a2010", linewidth=0.4),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color="#fc913a44", fill=NA, linewidth=0.5),
    axis.text        = element_text(color="#eae374", size=8),
    axis.title       = element_text(color="#fc913a", size=9),
    plot.title       = element_text(color="#f9d62e", size=11, face="bold"),
    legend.position  = "none",
    plot.margin      = margin(8,12,8,8))
}

# =============================================================================
# SERVER
# =============================================================================
server <- function(input, output, session) {
  
  use_metric      <- reactiveVal(TRUE)
  run_history     <- reactiveVal(list())
  results_store   <- reactiveVal(NULL)
  mc_store        <- reactiveVal(NULL)
  landing_pct_val <- reactiveVal(NULL)
  
  observeEvent(input$units, { use_metric(input$units == "metric") }, ignoreInit=TRUE)
  
  output$has_results <- reactive({ !is.null(results_store()) })
  output$has_mc      <- reactive({ !is.null(mc_store()) })
  outputOptions(output, "has_results", suspendWhenHidden=FALSE)
  outputOptions(output, "has_mc",      suspendWhenHidden=FALSE)
  
  # ---- SI stores (single source of truth; UI holds display units only) ----
  si <- list(
    dry_mass       = reactiveVal(0.090),
    diameter       = reactiveVal(0.024),
    body_length    = reactiveVal(0.300),
    cg_measured    = reactiveVal(0.220),
    nose_length    = reactiveVal(0.070),
    fin_root       = reactiveVal(0.050),
    fin_tip        = reactiveVal(0.025),
    fin_span       = reactiveVal(0.030),
    fin_sweep      = reactiveVal(0.020),
    fin_pos        = reactiveVal(0.320),
    parachute_diam = reactiveVal(0.305),
    rail_length    = reactiveVal(0.900),
    wind_speed     = reactiveVal(3.000),
    lug_length     = reactiveVal(0.025),
    lug_od         = reactiveVal(0.0040),
    lug_id         = reactiveVal(0.0032)
  )
  
  make_obs <- function(id, store, to_si, from_si) {
    prev_unit <- reactiveVal(NULL)
    observeEvent(input[[paste0(id,"_unit")]], {
      pu <- prev_unit(); nu <- input[[paste0(id,"_unit")]]
      if (!is.null(pu) && isTruthy(input[[id]])) {
        si_val <- to_si(input[[id]], pu); store(si_val)
        updateNumericInput(session, id, value=round(from_si(si_val, nu), 4))
      }
      prev_unit(nu)
    }, ignoreInit=FALSE)
    observeEvent(input[[id]], {
      u <- input[[paste0(id,"_unit")]]
      if (isTruthy(u) && isTruthy(input[[id]])) store(to_si(input[[id]], u))
    }, ignoreInit=TRUE)
  }
  make_length_observer <- function(id, store) make_obs(id, store, to_meters, from_meters)
  make_mass_observer   <- function(id, store) make_obs(id, store, to_kg,     from_kg)
  make_speed_observer  <- function(id, store) make_obs(id, store, to_ms,     from_ms)
  
  make_mass_observer(  "dry_mass_val",       si$dry_mass)
  make_length_observer("diameter",           si$diameter)
  make_length_observer("body_length",        si$body_length)
  make_length_observer("cg_measured",        si$cg_measured)
  make_length_observer("nose_length",        si$nose_length)
  make_length_observer("fin_root",           si$fin_root)
  make_length_observer("fin_tip",            si$fin_tip)
  make_length_observer("fin_span",           si$fin_span)
  make_length_observer("fin_sweep",          si$fin_sweep)
  make_length_observer("fin_pos",            si$fin_pos)
  make_length_observer("parachute_diameter", si$parachute_diam)
  make_length_observer("rail_length",        si$rail_length)
  make_speed_observer( "wind_speed_val",     si$wind_speed)
  make_length_observer("lug_length",         si$lug_length)
  make_length_observer("lug_od",             si$lug_od)
  make_length_observer("lug_id",             si$lug_id)
  
  # Write an SI value into a unit_input, respecting whatever unit is displayed
  set_si_length <- function(id, store, si_val) {
    if (!is.finite(si_val)) return(invisible(NULL))
    store(si_val)
    u <- input[[paste0(id,"_unit")]]; if (!isTruthy(u)) u <- "mm"
    updateNumericInput(session, id, value = round(from_meters(si_val, u), 4))
  }
  set_si_mass <- function(id, store, si_val) {
    if (!is.finite(si_val)) return(invisible(NULL))
    store(si_val)
    u <- input[[paste0(id,"_unit")]]; if (!isTruthy(u)) u <- "g"
    updateNumericInput(session, id, value = round(from_kg(si_val, u), 4))
  }
  
  motor_data <- reactive({
    if (!is.null(input$motor_file)) {
      tryCatch(parse_thrust_input(input$motor_file, NULL), error=function(e) NULL)
    } else {
      ec <- input$engine_choice
      if (is.null(ec) || ec == "") return(NULL)
      tryCatch(parse_thrust_input(NULL, ec), error=function(e) NULL)
    }
  })
  
  # ---- motor delay drives the ejection delay input ----
  observeEvent(motor_data(), {
    td <- motor_data()
    if (is.null(td) || length(td$delays) == 0) return()
    updateNumericInput(session, "parachute_delay", value = td$delays[1])
  })
  
  output$delay_hint <- renderUI({
    td <- motor_data()
    if (is.null(td)) return(tags$span(style="color:var(--dim);font-size:0.7rem;",
                                      "Select an engine to read its delay options."))
    if (length(td$delays) == 0)
      return(tags$span(style="color:var(--dim);font-size:0.7rem;",
                       "This motor is plugged / has no listed delay. Enter your own."))
    tags$span(style="color:var(--dim);font-size:0.7rem;",
              sprintf("Delays available on this motor: %s s. The charge fires this long AFTER burnout.",
                      paste(td$delays, collapse = ", ")))
  })
  
  aero_reactive <- reactive({
    if (!all(sapply(list(input$nose_type, input$fin_count), isTruthy))) return(NULL)
    if (!all(sapply(list(si$nose_length(), si$body_length(), si$diameter(),
                         si$fin_root(), si$fin_span(), si$cg_measured()),
                    function(x) isTruthy(x) && x > 0))) return(NULL)
    if (is.null(si$fin_tip()) || si$fin_tip() < 0) return(NULL)
    if ((si$fin_root() + si$fin_tip()) <= 0) return(NULL)
    tryCatch(
      compute_aero(input$nose_type,
                   si$nose_length(), si$body_length(), si$diameter(), input$fin_count,
                   si$fin_root(), si$fin_tip(), si$fin_span(), si$fin_sweep(),
                   si$fin_pos(), si$cg_measured(),
                   si$lug_length(), si$lug_od(), si$lug_id()),
      error=function(e) NULL)
  })
  
  # loaded / burnout stability
  aero_loaded <- reactive({
    aero <- aero_reactive(); if (is.null(aero)) return(NULL)
    td   <- motor_data()
    if (is.null(td)) {
      return(c(aero, list(stability_margin_loaded  = aero$stability_margin,
                          stability_margin_burnout = aero$stability_margin,
                          cg_loaded                = si$cg_measured())))
    }
    cg_motor     <- si$nose_length() + si$body_length() - td$motor_length_m / 2
    motor_mass   <- td$prop_mass + td$casing_mass
    total_loaded <- si$dry_mass() + motor_mass
    cg_loaded    <- (si$dry_mass() * si$cg_measured() + motor_mass * cg_motor) / total_loaded
    total_burnout <- si$dry_mass() + td$casing_mass
    cg_burnout    <- if (td$casing_mass > 0)
      (si$dry_mass() * si$cg_measured() + td$casing_mass * cg_motor) / total_burnout
    else
      si$cg_measured()
    c(aero, list(
      stability_margin_loaded  = (aero$CP - cg_loaded)  / si$diameter(),
      stability_margin_burnout = (aero$CP - cg_burnout) / si$diameter(),
      cg_loaded                = cg_loaded,
      cg_burnout               = cg_burnout
    ))
  })
  
  output$stability_indicator <- renderUI({
    al <- aero_loaded(); if (is.null(al)) return(NULL)
    sm  <- al$stability_margin_loaded
    fmt <- function(m) sprintf("%.1f mm  /  %.2f in", m*1000, m*39.3701)
    cls <- if (sm < 0.5) "stab-red" else if (sm < 1.0) "stab-yellow" else if (sm <= 3.0) "stab-green" else "stab-yellow"
    lbl <- if (sm < 0)   "UNSTABLE - CP AHEAD OF CG"
    else if (sm < 0.5) "UNSTABLE" else if (sm < 1.0) "MARGINAL"
    else if (sm <= 3.0) "STABLE"  else "OVERSTABLE"
    hint <- if (sm < 0.5) "Move CG forward or increase fin size. The simulator will fly this design as it really behaves." else
      if (sm > 3.0) "Risk of weathercocking in wind." else ""
    td <- motor_data()
    motor_lines <- if (!is.null(td)) {
      tagList(
        sprintf("CG loaded (ignition): %s", fmt(al$cg_loaded)), tags$br(),
        sprintf("Stability margin fully loaded:  %.2f cal", al$stability_margin_loaded), tags$br(),
        sprintf("Stability margin at burnout:    %.2f cal", al$stability_margin_burnout), tags$br(),
        tags$span(style="color:var(--dim);font-size:0.72rem;",
                  sprintf("(casing %.0f g, prop %.0f g)",
                          td$casing_mass*1000, td$prop_mass*1000))
      )
    } else {
      tags$span(style="color:var(--dim);", "Load an engine to see loaded CG")
    }
    div(class=paste("stab-box", cls),
        tags$b(sprintf("%s - %.2f cal (loaded)", lbl, sm)), tags$br(),
        sprintf("CP: %s", fmt(al$CP)), tags$br(),
        motor_lines,
        if (nchar(hint) > 0) tagList(tags$br(), tags$span(hint)) else NULL
    )
  })
  
  output$thrust_curve_plot <- renderPlot({
    tc <- motor_data()$thrust_curve
    if (is.null(tc) || nrow(tc) == 0) return(
      ggplot() + annotate("text", x=0.5, y=0.5, label="Select or upload an engine",
                          color="#6b7280", size=4) + theme_plot() +
        theme(axis.text=element_blank(), axis.title=element_blank(), panel.grid=element_blank()))
    tf    <- approxfun(tc$time, tc$thrust, yleft=0, yright=0)
    ti    <- integrate(tf, min(tc$time), max(tc$time))$value
    mt    <- max(tc$thrust)
    tc2   <- tc
    if (!use_metric()) tc2$thrust <- tc2$thrust * N_to_lbf
    ylab  <- if (use_metric()) "thrust (N)"       else "thrust (lbf)"
    t_ann <- if (use_metric()) sprintf("Total: %.2f Ns", ti) else sprintf("Total: %.2f lbf.s", ti*N_to_lbf)
    p_ann <- if (use_metric()) sprintf("Peak:  %.2f N",  mt) else sprintf("Peak:  %.2f lbf",   mt*N_to_lbf)
    ggplot(tc2, aes(time, thrust)) +
      geom_area(fill="#1a56db", alpha=0.08) +
      geom_line(color="#1a56db", linewidth=1) +
      geom_hline(yintercept=0, color="#e8eaf0") +
      annotate("text", x=max(tc2$time), y=max(tc2$thrust),
               label=t_ann, hjust=1, vjust=1.3, size=3.5, color="#f5f6f8") +
      annotate("text", x=max(tc2$time), y=max(tc2$thrust)*0.87,
               label=p_ann, hjust=1, vjust=1.3, size=3.5, color="#f5f6f8") +
      scale_y_continuous(limits=c(0, max(tc2$thrust)*1.18)) +
      labs(x="time (s)", y=ylab, title="Engine thrust vs. time") +
      theme_plot()
  })
  
  # ---- single simulation ----
  observeEvent(input$run, {
    aero   <- aero_reactive()
    parsed <- motor_data()
    if (is.null(aero))   { showNotification("Complete rocket geometry on Setup tab", type="warning"); return() }
    if (is.null(parsed)) { showNotification("Select an engine on Setup tab",         type="warning"); return() }
    sim <- tryCatch(
      flight_simulation_3d(
        parsed$thrust_curve, parsed$prop_mass, si$dry_mass(),
        parsed$casing_mass,
        si$diameter(), max(si$parachute_diam(), 0.05), input$parachute_delay,
        aero, aero$CNa_total, aero$CP, input$precision,
        si$wind_speed(), input$wind_dir,
        si$cg_measured(), si$nose_length(), si$body_length(),
        parsed$motor_length_m, si$rail_length(),
        input$launch_bearing, input$launch_angle,
        landing_only = FALSE,
        wind_turbulence_intensity = input$wind_turbulence_intensity,
        gust_duration = input$gust_duration),
      error=function(e) { showNotification(paste("Sim error:", e$message), type="error"); NULL })
    if (is.null(sim)) { showNotification("Simulation failed - check the .eng file", type="error"); return() }
    res  <- list(sim=sim, aero=aero,
                 burn_time = max(parsed$thrust_curve$time),
                 label=paste0("Run ", length(run_history())+1),
                 motor=if (isTruthy(input$engine_choice) && input$engine_choice!="") input$engine_choice else "custom")
    hist <- run_history(); hist[[length(hist)+1]] <- res; run_history(hist)
    results_store(res)
  })
  
  # ---- Monte Carlo ----
  observeEvent(input$run_mc, {
    aero   <- aero_reactive()
    parsed <- motor_data()
    if (is.null(aero))   { showNotification("Complete rocket geometry first", type="warning"); return() }
    if (is.null(parsed)) { showNotification("Select an engine first",         type="warning"); return() }
    n        <- input$mc_runs
    landings <- vector("list", n)
    withProgress(message="Monte Carlo", value=0, {
      for (i in seq_len(n)) {
        incProgress(1/n, detail=sprintf("Run %d / %d", i, n))
        sim <- tryCatch(
          flight_simulation_3d(
            parsed$thrust_curve,
            max(parsed$prop_mass * rnorm(1, 1, 0.01*input$prop_mass_std_dev), 1e-6),
            max(si$dry_mass()    * rnorm(1, 1, 0.01*input$dry_mass_std_dev),  1e-6),
            parsed$casing_mass,
            si$diameter(), max(si$parachute_diam(), 0.05),
            max(0, rnorm(1, input$parachute_delay, input$chute_delay_std_dev)),
            aero, aero$CNa_total, aero$CP,
            input$montecarlo_precision,
            max(0, rnorm(1, si$wind_speed(), 0.01*si$wind_speed()*input$wind_speed_std_dev)),
            rnorm(1, input$wind_dir, input$wind_dir_std_dev),
            si$cg_measured(), si$nose_length(), si$body_length(),
            parsed$motor_length_m, si$rail_length(),
            input$launch_bearing,
            max(0, rnorm(1, input$launch_angle, input$launch_angle_std_dev)),
            cd_scale     = max(rnorm(1, 1, 0.01 * input$cd_std_dev), 0.05),
            landing_only = TRUE,
            wind_turbulence_intensity = input$wind_turbulence_intensity,
            gust_duration = input$gust_duration),
          error=function(e) NULL)
        landings[[i]] <- if (!is.null(sim)) data.frame(x=sim$x, y=sim$y) else data.frame(x=0, y=0)
      }
    })
    mc_store(do.call(rbind, landings))
  })
  
  # ---- design save / load / .ork import ----
  design_list <- reactive({
    list(
      nose_type = input$nose_type, fin_count = input$fin_count,
      dry_mass = si$dry_mass(), diameter = si$diameter(),
      body_length = si$body_length(), cg_measured = si$cg_measured(),
      nose_length = si$nose_length(), fin_root = si$fin_root(),
      fin_tip = si$fin_tip(), fin_span = si$fin_span(),
      fin_sweep = si$fin_sweep(), fin_pos = si$fin_pos(),
      parachute_diam = si$parachute_diam(), rail_length = si$rail_length(),
      lug_length = si$lug_length(), lug_od = si$lug_od(), lug_id = si$lug_id(),
      engine_choice = input$engine_choice, parachute_delay = input$parachute_delay,
      launch_angle = input$launch_angle, launch_bearing = input$launch_bearing,
      units_note = "all lengths in metres, masses in kilograms"
    )
  })
  
  output$save_design <- downloadHandler(
    filename = function() paste0("rrrocket_design_", format(Sys.Date(), "%Y%m%d"), ".json"),
    content  = function(file) jsonlite::write_json(design_list(), file, auto_unbox = TRUE, pretty = TRUE)
  )
  
  apply_design <- function(d) {
    if (!is.null(d$nose_type))   updateSelectInput(session, "nose_type", selected = d$nose_type)
    if (!is.null(d$fin_count))   updateNumericInput(session, "fin_count", value = d$fin_count)
    if (!is.null(d$dry_mass))    set_si_mass("dry_mass_val", si$dry_mass, d$dry_mass)
    for (nm in c("diameter","body_length","cg_measured","nose_length","fin_root",
                 "fin_tip","fin_span","fin_sweep","fin_pos","rail_length",
                 "lug_length","lug_od","lug_id")) {
      if (!is.null(d[[nm]])) set_si_length(nm, si[[nm]], as.numeric(d[[nm]]))
    }
    if (!is.null(d$parachute_diam))
      set_si_length("parachute_diameter", si$parachute_diam, as.numeric(d$parachute_diam))
    if (!is.null(d$engine_choice))   updateSelectInput(session, "engine_choice", selected = d$engine_choice)
    if (!is.null(d$parachute_delay)) updateNumericInput(session, "parachute_delay", value = d$parachute_delay)
    if (!is.null(d$launch_angle))    updateNumericInput(session, "launch_angle", value = d$launch_angle)
    if (!is.null(d$launch_bearing))  updateSliderInput(session, "launch_bearing", value = d$launch_bearing)
  }
  
  observeEvent(input$load_design, {
    d <- tryCatch(jsonlite::read_json(input$load_design$datapath, simplifyVector = TRUE),
                  error = function(e) NULL)
    if (is.null(d)) { showNotification("Could not read that design file", type="error"); return() }
    apply_design(d)
    showNotification("Design loaded", type="message")
  })
  
  observeEvent(input$ork_file, {
    d <- tryCatch(parse_ork(input$ork_file$datapath), error = function(e) NULL)
    if (is.null(d) || length(d) == 0) {
      showNotification("Could not read that .ork file (only trapezoidal fin sets are supported)", type="error")
      return()
    }
    apply_design(d)
    showNotification("OpenRocket geometry imported. Check dry mass and CG - they are not imported reliably.",
                     type="warning", duration = 10)
  })
  
  or_overlay <- reactive({
    if (is.null(input$or_csv)) return(NULL)
    parse_or_csv(input$or_csv$datapath)
  })
  
  # ---- outputs ----
  output$run_history_table <- renderUI({
    hist <- run_history(); if (length(hist) == 0) return(p("No runs yet."))
    sc <- if (use_metric()) 1 else m_to_ft
    u  <- if (use_metric()) "m" else "ft"
    rows <- lapply(rev(seq_along(hist)), function(i) {
      r <- hist[[i]]; s <- r$sim
      tags$tr(tags$td(r$label), tags$td(r$motor),
              tags$td(sprintf("%.0f %s", max(s$altitude)*sc, u)),
              tags$td(sprintf("%.1f s",  max(s$time))),
              tags$td(sprintf("%.1f deg", ascent_max_aoa(s))),
              tags$td(sprintf("%.2f cal", min(s$stability_margin))))
    })
    tags$table(class="run-table",
               tags$thead(tags$tr(tags$th("Run"), tags$th("Motor"),
                                  tags$th("Apogee"), tags$th("Time"),
                                  tags$th("Max AoA"), tags$th("Min stab."))),
               tags$tbody(rows))
  })
  
  output$summary <- renderPrint({
    res <- results_store(); req(!is.null(res))
    r <- res$sim; ae <- res$aero
    sc <- if (use_metric()) 1 else m_to_ft
    u  <- if (use_metric()) "m"   else "ft"
    us <- if (use_metric()) "m/s" else "ft/s"
    
    t_ap  <- attr(r, "t_apogee"); t_ej <- attr(r, "t_eject")
    dv    <- attr(r, "deploy_v"); dalt <- attr(r, "deploy_alt")
    railv <- attr(r, "rail_exit_ms")
    
    cat(sprintf("apogee             %d %s\n",   round(max(r$altitude)*sc), u))
    cat(sprintf("max velocity       %.1f %s\n", max(r$velocity)*sc, us))
    if (!is.null(railv) && !is.na(railv)) {
      flag <- if (railv < V_RAIL_WARN) " *** LOW" else ""
      cat(sprintf("rail exit speed    %.1f %s%s\n", railv*sc, us, flag))
    }
    cat(sprintf("max Mach           %.3f\n", max(r$mach)))
    cat(sprintf("max AoA (to apogee) %.1f deg\n", ascent_max_aoa(r)))
    cat(sprintf("time to apogee     %.2f s\n", r$time[which.max(r$altitude)]))
    cat(sprintf("total flight time  %.2f s\n", max(r$time)))
    cat("\n")
    cat(sprintf("burnout            %.2f s\n", attr(r, "burn_time")))
    cat(sprintf("ejection fires     %.2f s\n", t_ej))
    if (!is.na(t_ap)) {
      d <- t_ej - t_ap
      cat(sprintf("  vs apogee        %.2f s %s\n", abs(d),
                  if (d < -0.05) "EARLY" else if (d > 0.05) "LATE" else "(on apogee)"))
    }
    if (!is.na(dv)) {
      cat(sprintf("deployment speed   %.1f %s%s\n", dv*sc, us,
                  if (dv > V_DEPLOY_WARN) " *** HIGH" else ""))
      cat(sprintf("deployment alt     %d %s\n", round(dalt*sc), u))
    } else {
      cat("deployment         *** NEVER FIRED BEFORE IMPACT ***\n")
    }
    cat("\n")
    cat(sprintf("SM at ignition     %.2f cal\n", r$stability_margin[1]))
    cat(sprintf("SM minimum         %.2f cal\n", min(r$stability_margin)))
    cat(sprintf("Cd0 nominal        %.4f   (50 m/s, M=0)\n", ae$Cd))
    cat(sprintf("  nose             %.4f\n", ae$Cd_nose))
    cat(sprintf("  body             %.4f\n", ae$Cd_body))
    cat(sprintf("  fins             %.4f\n", ae$Cd_fins))
    cat(sprintf("  base             %.4f\n", ae$Cd_base))
    cat(sprintf("  launch lug       %.4f\n", ae$Cd_lug))
    
    ov <- or_overlay()
    if (!is.null(ov)) {
      cat("\n")
      cat(sprintf("OpenRocket apogee  %d %s\n", round(max(ov$altitude)*sc), u))
      cat(sprintf("  difference       %+.1f %%\n",
                  100 * (max(r$altitude) - max(ov$altitude)) / max(ov$altitude)))
    }
  })
  
  output$safety_box <- renderUI({
    res <- results_store(); req(!is.null(res)); r <- res$sim
    warns <- character(0)
    railv <- attr(r, "rail_exit_ms")
    dv    <- attr(r, "deploy_v")
    t_ap  <- attr(r, "t_apogee"); t_ej <- attr(r, "t_eject")
    
    if (is.na(dv))
      warns <- c(warns, sprintf("Ballistic impact warning. The ejection charge fires at %.1f s, which is after the rocket hits the ground at %.1f s.", t_ej, max(r$time)))
    if (!is.na(railv) && railv < V_RAIL_WARN)
      warns <- c(warns, sprintf("Rail exit speed is %.1f m/s. Below ~%d m/s the fins cannot fully stabilize the rocket and prevent weathercocking..", railv, V_RAIL_WARN))
    if (!is.na(dv) && dv > V_DEPLOY_WARN)
      warns <- c(warns, sprintf("Chute opens at %.1f m/s. High-speed deployment is potentially dangerous.", dv))
    if (!is.na(t_ap) && (t_ej - t_ap) < -0.5)
      warns <- c(warns, sprintf("Ejection fires %.1f s before apogee.", t_ap - t_ej))
    if (!is.na(t_ap) && (t_ej - t_ap) > 1.5)
      warns <- c(warns, sprintf("Ejection fires %.1f s after apogee.", t_ej - t_ap))
    if (min(r$stability_margin) < 1.0)
      warns <- c(warns, sprintf("Minimum in-flight stability margin %.2f cal. Below 1 cal the rocket is only marginally stable.", min(r$stability_margin)))
    if (ascent_max_aoa(r) > 30)
      warns <- c(warns, sprintf("Max angle of attack %.0f deg during the ascent.", ascent_max_aoa(r)))
    if (ascent_max_aoa(r) > 12 && ascent_max_aoa(r) <= 30)
      warns <- c(warns, sprintf("Max angle of attack %.0f deg during the ascent.", ascent_max_aoa(r)))
    
    if (length(warns) == 0)
      return(div(class="stab-box stab-green", tags$b("No safety flags on this flight.")))
    div(class="stab-box stab-red",
        tags$b("Safety flags"), tags$br(),
        tagList(lapply(warns, function(w) tagList(tags$span(paste0("- ", w)), tags$br()))))
  })
  
  output$fin_preview <- renderPlot({
    root  <- si$fin_root(); tip <- si$fin_tip(); span <- si$fin_span()
    sweep <- si$fin_sweep(); diam <- si$diameter()
    
    req(isTruthy(root) && root > 0,
        isTruthy(tip)  && tip  >= 0,
        isTruthy(span) && span > 0,
        isTruthy(sweep)&& sweep >= 0,
        isTruthy(diam) && diam > 0)
    
    du  <- if (isTruthy(input$fin_root_unit)) input$fin_root_unit else "mm"
    fmt <- function(v) sprintf("%.1f %s", from_meters(v, du), du)
    
    body_r <- diam / 2
    bx     <- -body_r
    
    rx1 <- 0;    ry1 <- 0
    rx2 <- 0;    ry2 <- root
    tx1 <- span; ty1 <- sweep
    tx2 <- span; ty2 <- sweep + tip
    
    fin_df <- data.frame(x = c(rx1, tx1, tx2, rx2), y = c(ry1, ty1, ty2, ry2))
    
    pad_x    <- span * 0.55
    pad_y    <- max(root, sweep + tip) * 0.32
    xlim     <- c(bx - body_r * 0.3, span + pad_x)
    ylim     <- c(-pad_y, max(root, sweep + tip) + pad_y)
    
    off_root <- -span * 0.1
    off_tip  <-  span * 0.05
    off_span <-  max(root, sweep + tip) * 0.18
    
    ann_col  <- "#fff8f0"; fin_fill <- "#eae37433"; fin_col <- "#f9d62e"
    body_col <- "#3a2510"; dim_col  <- "#1a56db";   txt_col <- "#fff8f0"
    
    ggplot() +
      annotate("rect", xmin = bx, xmax = 0,
               ymin = -pad_y * 0.6, ymax = max(root, sweep + tip) + pad_y * 0.6,
               fill = body_col, color = "#eae374aa", linewidth = 0.6) +
      geom_polygon(data = fin_df, aes(x = x, y = y),
                   fill = fin_fill, color = fin_col, linewidth = 1.2) +
      annotate("segment", x = off_root, xend = off_root, y = ry1, yend = ry2,
               color = dim_col, linewidth = 0.7,
               arrow = arrow(ends = "both", length = unit(5, "pt"), type = "closed")) +
      annotate("segment", x = off_root - span*0.03, xend = off_root + span*0.01,
               y = ry1, yend = ry1, color = dim_col, linewidth = 0.5) +
      annotate("segment", x = off_root - span*0.03, xend = off_root + span*0.01,
               y = ry2, yend = ry2, color = dim_col, linewidth = 0.5) +
      annotate("text", x = off_root - span * 0.06, y = (ry1 + ry2) / 2,
               label = paste0("root\n", fmt(root)),
               color = txt_col, size = 3.1, hjust = 1, fontface = "bold") +
      annotate("segment", x = span + off_tip, xend = span + off_tip, y = ty1, yend = ty2,
               color = dim_col, linewidth = 0.7,
               arrow = arrow(ends = "both", length = unit(5, "pt"), type = "closed")) +
      annotate("segment", x = span + off_tip - span*0.01, xend = span + off_tip + span*0.04,
               y = ty1, yend = ty1, color = dim_col, linewidth = 0.5) +
      annotate("segment", x = span + off_tip - span*0.01, xend = span + off_tip + span*0.04,
               y = ty2, yend = ty2, color = dim_col, linewidth = 0.5) +
      annotate("text", x = span + off_tip + span * 0.06, y = (ty1 + ty2) / 2,
               label = paste0("tip\n", fmt(tip)),
               color = txt_col, size = 3.1, hjust = 0, fontface = "bold") +
      annotate("segment", x = 0, xend = span, y = -off_span, yend = -off_span,
               color = dim_col, linewidth = 0.7,
               arrow = arrow(ends = "both", length = unit(5, "pt"), type = "closed")) +
      annotate("text", x = span / 2, y = -off_span - 0.0024,
               label = paste0("semi-span  ", fmt(span)),
               color = txt_col, size = 3.1, hjust = 0.5, fontface = "bold") +
      { if (sweep > 1e-5) list(
        annotate("segment", x = span + off_tip, xend = span + off_tip, y = ry1, yend = ty1,
                 color = dim_col, linewidth = 0.6,
                 arrow = arrow(ends = "both", length = unit(4, "pt"), type = "closed")),
        annotate("segment", x = rx1, xend = tx1, y = ty1, yend = ty1,
                 color = dim_col, linewidth = 0.4, linetype = "dotted"),
        annotate("segment", x = tx1, xend = tx1, y = ry1, yend = ty1,
                 color = dim_col, linewidth = 0.4, linetype = "dotted"),
        annotate("text", x = span * 1.13, y = sweep * 0.3,
                 label = paste0("sweep\n", fmt(sweep)),
                 color = ann_col, size = 2.9, hjust = 0, fontface = "bold")
      ) else list() } +
      annotate("text", x = bx / 2, y = max(root, sweep + tip) + pad_y * 0.7,
               label = "body wall", color = "#eae374aa", size = 2.6, hjust = 0.5) +
      scale_x_continuous(expand = expansion(0)) +
      scale_y_continuous(expand = expansion(0)) +
      coord_fixed(xlim = xlim, ylim = ylim) +
      labs(title = "Fin preview", x = NULL, y = NULL) +
      theme_plot() +
      theme(axis.text = element_blank(), axis.ticks = element_blank(),
            panel.grid = element_blank())
  }, bg = "#120800")
  
  output$altitude_plot <- renderPlot({
    res <- results_store(); req(!is.null(res)); r <- res$sim
    sc   <- if (use_metric()) 1 else m_to_ft
    ylab <- if (use_metric()) "altitude (m)" else "altitude (ft)"
    p <- ggplot(data.frame(t=r$time, alt=r$altitude*sc), aes(t, alt)) +
      geom_area(fill="#1a56db", alpha=0.15) +
      geom_line(color="#1a56db", linewidth=1) +
      geom_hline(yintercept=0, color="#e8eaf0")
    ov <- or_overlay()
    if (!is.null(ov)) {
      p <- p + geom_line(data = data.frame(t = ov$time, alt = ov$altitude*sc),
                         aes(t, alt), color="#ff4e50", linewidth=0.8, linetype="dashed")
    }
    ttl <- if (is.null(ov)) "Altitude" else "Altitude (dashed red = OpenRocket)"
    p + labs(x="time (s)", y=ylab, title=ttl) + theme_plot()
  })
  
  output$velocity_plot <- renderPlot({
    res <- results_store(); req(!is.null(res)); r <- res$sim
    sc   <- if (use_metric()) 1 else m_to_ft
    ylab <- if (use_metric()) "vertical velocity (m/s)" else "vertical velocity (ft/s)"
    ggplot(data.frame(t=r$time, vz=r$vz*sc), aes(t, vz)) +
      geom_line(color="#f9d62e", linewidth=1) +
      geom_hline(yintercept=0, color="#fc913a44", linetype="dashed") +
      labs(x="time (s)", y=ylab, title="Vertical velocity") + theme_plot()
  })
  
  output$aoa_plot <- renderPlot({
    res <- results_store(); req(!is.null(res)); r <- res$sim
    d <- r[ascent_rows(r), ]
    ggplot(d, aes(time, aoa)) +
      geom_line(color="#ff4e50", linewidth=1) +
      geom_hline(yintercept=17, color="#fc913a66", linetype="dashed") +
      labs(x="time (s)", y="angle of attack (deg)",
           title="Angle of attack (liftoff to apogee)") + theme_plot()
  })
  
  output$stab_plot <- renderPlot({
    res <- results_store(); req(!is.null(res)); r <- res$sim
    d <- r[ascent_rows(r), ]
    ggplot(d, aes(time, stability_margin)) +
      geom_line(color="#e2f4c7", linewidth=1) +
      geom_hline(yintercept=1, color="#22c55e88", linetype="dashed") +
      geom_hline(yintercept=0, color="#ef444488") +
      labs(x="time (s)", y="stability margin (cal)",
           title="Live stability margin") + theme_plot()
  })
  
  output$track_3d <- renderPlotly({
    res <- results_store(); req(!is.null(res)); r <- res$sim
    sc <- if (use_metric()) 1 else m_to_ft
    xl <- if (use_metric()) "East (m)"     else "East (ft)"
    yl <- if (use_metric()) "North (m)"    else "North (ft)"
    zl <- if (use_metric()) "Altitude (m)" else "Altitude (ft)"
    
    phase_colors <- c("1" = "#ff4e50", "2" = "#f9d62e", "3" = "#e2f4c7")
    phase_names  <- c("1" = "Boost",   "2" = "Coast",   "3" = "Descent")
    
    r$segment <- cumsum(c(1, diff(r$phase) != 0))
    
    traces <- lapply(unique(r$segment), function(seg) {
      d  <- r[r$segment == seg, ]
      ph <- as.character(d$phase[1])
      next_row <- r[r$segment == seg + 1, ]
      if (nrow(next_row) > 0) d <- rbind(d, next_row[1, ])
      list(x = d$x * sc, y = d$y * sc, z = d$altitude * sc,
           col = phase_colors[ph], name = phase_names[ph], ph = ph)
    })
    
    fig <- plot_ly(type = "scatter3d", mode = "lines")
    seen_phases <- character(0)
    for (tr in traces) {
      show_legend <- !(tr$ph %in% seen_phases)
      seen_phases <- union(seen_phases, tr$ph)
      fig <- fig |> add_trace(
        x = tr$x, y = tr$y, z = tr$z,
        type = "scatter3d", mode = "lines",
        name = tr$name, showlegend = show_legend,
        line = list(color = tr$col, width = 4))
    }
    
    fig <- fig |>
      add_trace(
        x = tail(r$x, 1) * sc, y = tail(r$y, 1) * sc, z = 0,
        type = "scatter3d", mode = "markers", name = "Landing",
        marker = list(color = "#e02424", size = 7, symbol = "circle", opacity = 0),
        projection = list(z = list(show = TRUE, opacity = 1, scale = 1)),
        showlegend = TRUE) |>
      add_trace(
        x = tail(r$x, 1) * sc, y = tail(r$y, 1) * sc, z = 0,
        type = "scatter3d", mode = "markers", name = "Landing",
        marker = list(color = "rgba(0,0,0,0)", size = 14, symbol = "circle-open",
                      opacity = 0, line = list(color = "#ffffff", width = 2)),
        projection = list(z = list(show = TRUE, opacity = 1, scale = 1)),
        showlegend = FALSE) |>
      add_trace(
        x = tail(r$x, 1) * sc, y = tail(r$y, 1) * sc, z = 0,
        type = "scatter3d", mode = "markers", name = "Landing",
        marker = list(color = "rgba(0,0,0,0)", size = 5, symbol = "circle",
                      opacity = 0, line = list(color = "#e02424", width = 0)),
        projection = list(z = list(show = TRUE, opacity = 1, scale = 1)),
        showlegend = FALSE) |>
      layout(
        paper_bgcolor = "#120800",
        font  = list(color = "#f9d62e", family = "Lexend, sans-serif"),
        legend = list(bgcolor = "rgba(26,16,8,0.85)", bordercolor = "#fc913a44",
                      borderwidth = 1, font = list(color = "#eae374", size = 11)),
        scene = list(
          bgcolor = "#1a0a02",
          xaxis = list(title = xl, gridcolor = "#3a2010", color = "#eae374"),
          yaxis = list(title = yl, gridcolor = "#3a2010", color = "#eae374"),
          zaxis = list(title = zl, gridcolor = "#3a2010", color = "#eae374")))
    fig
  })
  
  # ---- map / Monte Carlo landing zone ----
  launch_point  <- reactiveVal(list(lat=38.89, lng=-77.03))
  drawn_polygon <- reactiveVal(NULL)
  
  observeEvent(input$map_click,            { launch_point(list(lat=input$map_click$lat, lng=input$map_click$lng)) })
  observeEvent(input$map_draw_new_feature, { drawn_polygon(input$map_draw_new_feature) })
  
  rocket_icon <- makeIcon(
    iconUrl="https://raw.githubusercontent.com/pointhi/leaflet-color-markers/master/img/marker-icon-2x-red.png",
    iconWidth=25, iconHeight=41, iconAnchorX=12, iconAnchorY=41)
  
  output$map <- renderLeaflet({
    leaflet() |> addProviderTiles("Esri.WorldImagery") |>
      setView(lat=36.9052, lng=-81.0768, zoom=5) |>
      addDrawToolbar(polylineOptions=FALSE, circleOptions=drawCircleOptions(),
                     markerOptions=FALSE, circleMarkerOptions=FALSE,
                     rectangleOptions=drawRectangleOptions(),
                     polygonOptions=drawPolygonOptions(), editOptions=FALSE)
  })
  
  observe({
    lp <- launch_point()
    leafletProxy("map") |> clearGroup("launch") |>
      addMarkers(lng=lp$lng, lat=lp$lat, icon=rocket_icon, label="Launch pad", group="launch")
  })
  
  output$landing_pct <- renderPrint({
    p <- landing_pct_val(); req(!is.null(p))
    lc <- mc_store()
    cat(sprintf("%.1f%% of %d simulated flights land inside the polygon\n", p, nrow(lc)))
    cat(sprintf("mean landing distance from pad: %.0f m\n",
                mean(sqrt(lc$x^2 + lc$y^2))))
    cat(sprintf("95th percentile distance:       %.0f m\n",
                quantile(sqrt(lc$x^2 + lc$y^2), 0.95)))
  })
  
  observeEvent(list(mc_store(), drawn_polygon()), {
    lc <- mc_store(); poly <- drawn_polygon()
    if (is.null(lc) || is.null(poly)) return()
    lp   <- launch_point(); lat0 <- lp$lat; lng0 <- lp$lng
    lc$lat <- lat0 + (lc$y / 111320)
    lc$lng <- lng0 + (lc$x / (111320 * cos(lat0 * pi/180)))
    poly_coords <- poly$geometry$coordinates[[1]]
    poly_mat    <- do.call(rbind, lapply(poly_coords, function(p) c(p[[1]], p[[2]])))
    if (!identical(poly_mat[1,], poly_mat[nrow(poly_mat),]))
      poly_mat <- rbind(poly_mat, poly_mat[1,])
    pip <- function(px, py, pm) sapply(seq_along(px), function(k) {
      x<-px[k]; y<-py[k]; n<-nrow(pm); j<-n; inside<-FALSE
      for (i in 1:n) {
        xi<-pm[i,1]; yi<-pm[i,2]; xj<-pm[j,1]; yj<-pm[j,2]
        if (((yi>y) != (yj>y)) && (x < (xj-xi)*(y-yi)/(yj-yi)+xi))
          inside <- !inside
        j <- i
      }
      inside
    })
    inside <- pip(lc$lng, lc$lat, poly_mat)
    pct    <- mean(inside) * 100
    landing_pct_val(pct)
    lc$color <- ifelse(inside, "#0e9f6e", "#e02424")
    leafletProxy("map") |> clearMarkers() |> clearGroup("launch") |>
      addCircleMarkers(lng=lc$lng, lat=lc$lat, radius=4, color=lc$color,
                       fillOpacity=0.7, stroke=FALSE) |>
      addMarkers(lng=lng0, lat=lat0, icon=rocket_icon, label="Launch pad", group="launch") |>
      addPopups(lng=lng0, lat=lat0, popup=sprintf("%.1f%% land inside safe zone", pct))
  })
}

shinyApp(ui=ui, server=server)