library(shiny)
library(bslib)
library(plotly)
library(ggplot2)
library(leaflet)
library(leaflet.extras)

unit_choices_length <- c("mm", "cm", "in", "ft", "m")
unit_choices_mass   <- c("g", "oz", "kg", "lb")
unit_choices_speed  <- c("m/s", "mph", "km/h", "knots")

# modeled off of OpenRocket technical documentation
# Nisanken 2013, Barrowman 1967

g0      <- 9.80665
R_air   <- 287.058
gamma_a <- 1.4
Cd_chute <- 0.80 # from OpenRocket
m_to_ft  <- 3.28084
N_to_lbf <- 0.224809

#unit conversions
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

# only troposphere
isa_atm <- function(z) {
  z   <- max(z, 0)
  T   <- 288.15 - 0.0065 * z
  rho <- 1.225 * (T / 288.15)^4.2561
  list(rho = rho, a = sqrt(gamma_a * R_air * T))
}

# implement Cf calcs from openrocket
# Turbulent smooth:     Cf = 1 / (1.50*ln(R) - 5.6)^2
# Roughness-limited:    Cf = 0.032 * (Rs/L)^0.2
# Low-Re floor:         Cf = 0.0148  (R < 1e4)
# Subsonic Mach corr.:  Cf_c = Cf * (1 - 0.1*M^2)        [eq 3.82]
# Supersonic Mach corr: Cf_c = Cf / (1 + 0.15*M^2)^0.58  [eq 3.83]
# Rs = surface roughness height (m).  Default 60e-6 m ("optimum
#      paint-sprayed" finish, mid-range of Table 3.2) — the caller
#      can override.

skin_friction_cf <- function(velocity, char_length, mach,
                             Rs = 60e-6) {
  # Kinematic viscosity of air at ~15 °C (ISA sea level)
  nu <- 1.461e-5
  Re <- max(velocity * char_length / nu, 1)
  
  # Critical Reynolds number for roughness transition (eq 3.79):
  #   Rcrit = 51 * (Rs/L)^(-1.039)
  
  Rcrit <- 51 * (Rs / char_length)^(-1.039)
  
  if (Re < 1e4) {
    Cf <- 1.48e-2                              # low-Re floor
  } else if (Re < Rcrit) {
    Cf <- 1 / (1.50 * log(Re) - 5.6)^2        # turbulent smooth (eq 3.78)
  } else {
    Cf <- 0.032 * (Rs / char_length)^0.2      # roughness-limited (eq 3.80)
  }
  
  # air compresses around rocket, reduces skin friction at higher mach#
  if (mach < 1.0) {
    Cf_c <- Cf * (1.0 - 0.1 * mach^2)
  } else {
    Cf_c <- Cf / (1.0 + 0.15 * mach^2)^0.58
    if (Re >= Rcrit) {
      Cf_c_rough <- Cf / (1.0 + 0.18 * mach^2)
      Cf_c <- max(Cf_c, Cf_c_rough)
    }
  }
  
  Cf_c
}

# multiply Cd by these numbers for each mach #
#   M < 0.8 : Prandtl-Glauert  1/sqrt(1-M^2)
#   0.8–1.0 : continuous transonic ramp → 2.4 at M=1
#   1.0–2.0 : supersonic decay → 1.6 at M=2
#   M ≥ 2.0 : ~constant 1.6
mach_cd_factor <- function(M) {
  if      (M < 0.8) 1 / sqrt(max(1 - M^2, 0.01))
  else if (M < 1.0) 1.6667 + 3.6667 * (M - 0.8)
  else if (M < 2.0) 2.4 - 0.8 * (M - 1.0)
  else              1.6
}

# Openrocket
compute_aero <- function(nose_type, nose_length, body_length,
                         bt_diameter, fin_count, fin_root,
                         fin_tip, fin_span, fin_sweep, cg_measured,
                         ref_velocity = 50) {
  
  bt_radius <- bt_diameter / 2
  Aref      <- pi * bt_radius^2
  
  Xcp_nose <- switch(nose_type,
                     conical    = (2/3) * nose_length,
                     ogive      = 0.466 * nose_length,
                     parabolic  = 0.5   * nose_length)
  
  CNa_nose <- 2.0
  
  ha <- atan(bt_radius / nose_length)
  Cd_nose_pressure <- switch(nose_type,
                             conical   = 0.8 * sin(ha)^2,
                             ogive     = 0.5 * sin(ha)^2,
                             parabolic = 0.3 * (bt_diameter / nose_length)^2)
  
  nose_slant <- sqrt(nose_length^2 + bt_radius^2)
  Awet_nose  <- pi * bt_radius * nose_slant
  rocket_length <- nose_length + body_length
# nose
  Cf_nose <- skin_friction_cf(ref_velocity, rocket_length, mach = 0)
  Cd_nose_friction <- Cf_nose * Awet_nose / Aref
  Cd_nose <- Cd_nose_pressure + Cd_nose_friction
# bt
  Awet_body <- pi * bt_diameter * body_length
  fB        <- rocket_length / bt_diameter          # fineness ratio
  Cf_body   <- skin_friction_cf(ref_velocity, rocket_length, mach = 0)
  Cd_body   <- Cf_body * (1 + 2 / fB) * Awet_body / Aref
  
# fins
  Afin_one  <- 0.5 * (fin_root + fin_tip) * fin_span
  Kfb       <- 1 + bt_radius / (fin_span + bt_radius)
  CNa_fin   <- Kfb * (4 * fin_count * (fin_span / bt_diameter)^2) /
    (1 + sqrt(1 + (2 * fin_span / (fin_root + fin_tip))^2))
  
  Xb      <- nose_length + body_length
  Xcp_fin <- Xb +
    (fin_sweep / 3) * ((fin_root + 2 * fin_tip) / (fin_root + fin_tip)) +
    (1 / 6) * (fin_root + fin_tip - fin_root * fin_tip / (fin_root + fin_tip))
  
  Awet_fins <- 2 * fin_count * Afin_one
 
  # calc fin area
  c_bar_fin <- (fin_root + fin_tip) / 2
  Cf_fins   <- skin_friction_cf(ref_velocity, c_bar_fin, mach = 0)
  Cd_fins <- Cf_fins * Awet_fins / Aref
  
  # find center of pressure
  CNa_total <- CNa_nose + CNa_fin
  CP        <- (CNa_nose * Xcp_nose + CNa_fin * Xcp_fin) / CNa_total
  
  # cd_base only used to calculate stability for popup, not for sim
  Cd_base <- 0.12   # M = 0 reference value
  
  Cd_parasite <- Cd_nose + Cd_body + Cd_fins
  Cd_total    <- Cd_parasite + Cd_base
  
  list(
    Cd               = Cd_total,
    Cd_nose          = Cd_nose,
    Cd_body          = Cd_body,
    Cd_fins          = Cd_fins,
    Cd_base          = Cd_base,
    Cd_parasite      = Cd_parasite,
    CP               = CP,
    CNa_total        = CNa_total,
    stability_margin = (CP - cg_measured) / bt_diameter
  )
}

# cd during flight

live_base_drag <- function(M) {
  if (M < 1.0) 0.12 + 0.13 * M^2
  else         0.25 / M
}

# approximates 6DOF weathercocking, saves computational complexity
weathercock_gain <- function(sm) max(min(sm / 3.0, 1.0) * 0.6, 0)

###sim
flight_simulation_3d <- function(thrust_curve, prop_mass, dry_mass,
                                 casing_mass,
                                 bt_diameter, chute_diameter, chute_delay,
                                 Cd_parasite, CNa_total, CP, precision,
                                 wind_speed_ref, wind_dir_deg,
                                 cg_dry_m, nose_length, body_length,
                                 motor_length_m, rail_length_m,
                                 launch_bearing_deg, launch_angle_deg,
                                 landing_only = FALSE) {
  
  if (is.null(thrust_curve) || nrow(thrust_curve) == 0) return(NULL)
  
  thrust    <- approxfun(thrust_curve$time, thrust_curve$thrust,
                         yleft = 0, yright = 0)
  burn_time <- max(thrust_curve$time)
  
  eff_dry_mass <- dry_mass + casing_mass
  
  total_imp <- integrate(thrust,
                         min(thrust_curve$time),
                         max(thrust_curve$time))$value
  v_exhaust <- total_imp / prop_mass
  if (!is.finite(v_exhaust) || v_exhaust <= 0) {
    showNotification("Invalid .eng file: check thrust curve", type = "error")
    return(NULL)
  }
  
  bt_area    <- pi * (bt_diameter / 2)^2
  chute_area <- pi * (chute_diameter / 2)^2
  
  cg_motor <- nose_length + body_length - motor_length_m / 2
  
  ar <- launch_angle_deg   * pi / 180
  br <- launch_bearing_deg * pi / 180
  lx <- sin(ar) * sin(br)
  ly <- sin(ar) * cos(br)
  lz <- cos(ar)
  
  rail_ht <- rail_length_m * cos(ar)

  # Euler's method to keep runs short and feasible for mc
  x  <- 0; y  <- 0; z  <- 0
  vx <- 0; vy <- 0; vz <- 0
  m  <- eff_dry_mass + prop_mass
  m_prop <- prop_mass
  t  <- 0
  t_apogee   <- NA
  chute_open <- FALSE
  on_rail    <- TRUE
  
  wd_rad <- wind_dir_deg * pi / 180
  wx <- -wind_speed_ref * sin(wd_rad)
  wy <- -wind_speed_ref * cos(wd_rad)
  
  max_steps <- ceiling(1200 / max(precision, 0.001))
  if (!landing_only) {
    out <- data.frame(
      time             = numeric(max_steps),
      x                = numeric(max_steps),
      y                = numeric(max_steps),
      altitude         = numeric(max_steps),
      velocity         = numeric(max_steps),
      vx               = numeric(max_steps),
      vy               = numeric(max_steps),
      vz               = numeric(max_steps),
      mach             = numeric(max_steps),
      stability_margin = numeric(max_steps)
    )
    i <- 1L
  }
  
  repeat {
    dt <- if (t <= burn_time + 1.0) min(precision, 0.01) else max(precision, 0.02)
    
    atm   <- isa_atm(z)
    rho   <- atm$rho
    a_snd <- atm$a
    
# wind physics, gusts follow Normal dist, 
    z_ref <- 10.0
    shear <- (max(z, z_ref) / z_ref)^0.14
    mu_x  <- -wind_speed_ref * shear * sin(wd_rad)
    mu_y  <- -wind_speed_ref * shear * cos(wd_rad)
    

    local_wind_speed <- wind_speed_ref * shear
    sigma_w <- 0.15 * max(local_wind_speed, 0.01)
    
    # OU update: mean-revert toward altitude-corrected mean, add Gaussian gust.
    # Decorrelation rate alpha = 0.5 s^-1
    wx <- wx + 0.5 * (mu_x - wx) * dt + sigma_w * sqrt(dt) * rnorm(1)
    wy <- wy + 0.5 * (mu_y - wy) * dt + sigma_w * sqrt(dt) * rnorm(1)
    
    # relative vel in each coordinate
    vrx <- vx - wx
    vry <- vy - wy
    vrz <- vz
    vrm <- max(sqrt(vrx^2 + vry^2 + vrz^2), 1e-6)
    M   <- vrm / a_snd
    
    # drag force
    if (chute_open) {
      Cd_eff   <- Cd_chute   # [OR-FIX] 0.80 (Hoerner default, matches OpenRocket)
      area_eff <- chute_area
    } else {
      Cd_eff   <- Cd_parasite * mach_cd_factor(M) + live_base_drag(M)
      area_eff <- bt_area
    }
    q  <- 0.5 * rho * vrm^2
    Fd <- Cd_eff * area_eff * q
    Fdx <- -Fd * (vrx / vrm)
    Fdy <- -Fd * (vry / vrm)
    Fdz <- -Fd * (vrz / vrm)
    
    Ft <- thrust(t)
    
    cg_live <- (eff_dry_mass * cg_dry_m + m_prop * cg_motor) / m
    sm_live <- (CP - cg_live) / bt_diameter
    
    # weathercocking side force
    Fwx <- 0; Fwy <- 0
    if (on_rail && z >= rail_ht) on_rail <- FALSE
    if (!on_rail && is.na(t_apogee) && vz > 0.5) {
      lat_wx  <- wx - vx
      lat_wy  <- wy - vy
      lat_wm  <- max(sqrt(lat_wx^2 + lat_wy^2), 1e-6)
      sin_aoa <- min(lat_wm / vrm, 1.0)
      wc      <- weathercock_gain(sm_live)
      Fw      <- wc * CNa_total * q * bt_area * sin_aoa
      Fwx     <- Fw * (lat_wx / lat_wm)
      Fwy     <- Fw * (lat_wy / lat_wm)
    }
    
    # kinematics
    ax <- (Fdx + Fwx + Ft * lx) / m
    ay <- (Fdy + Fwy + Ft * ly) / m
    az <- (Ft  * lz  + Fdz - m * g0) / m
    
    if (!landing_only) {
      out$time[i]             <- t
      out$x[i]                <- x
      out$y[i]                <- y
      out$altitude[i]         <- z
      out$velocity[i]         <- sqrt(vx^2 + vy^2 + vz^2)
      out$vx[i]               <- vx
      out$vy[i]               <- vy
      out$vz[i]               <- vz
      out$mach[i]             <- M
      out$stability_margin[i] <- sm_live
      i <- i + 1L
    }
    
    if (t <= burn_time) {
      m_prop <- max(m_prop - (Ft / v_exhaust) * dt, 0)
      m      <- eff_dry_mass + m_prop
    }
    if (!is.na(t_apogee) && !chute_open && (t - t_apogee) >= chute_delay)
      chute_open <- TRUE
    
    vx <- vx + ax * dt;  vy <- vy + ay * dt;  vz <- vz + az * dt
    x  <- x  + vx * dt;  y  <- y  + vy * dt;  z  <- z  + vz * dt
    t  <- t  + dt
    
    if (is.na(t_apogee) && z > 1 && vz < 0) t_apogee <- t
    
    if (!is.na(t_apogee) && z <= 0) {
      if (landing_only) {
        frac <- if (abs(vz * dt) > 1e-9) z / (vz * dt) else 0
        x <- x - vx * dt * frac
        y <- y - vy * dt * frac
        return(data.frame(x = x, y = y))
      }
      break
    }
    if (t > 1200 || (!landing_only && i >= max_steps)) break
  }
  
  if (landing_only) return(data.frame(x = x, y = y))
  
  result <- out[1:(i - 1), ]
  rail_row <- which(result$altitude >= rail_ht)
  attr(result, "rail_exit_ms") <-
    if (length(rail_row) > 0) result$velocity[rail_row[1]] else NA_real_
  result
}

# ── Parse .eng motor file ────────────────────────────────────────────────────
parse_thrust_input <- function(motor_file, engine_choice) {
  read_eng <- function(lines) {
    lines <- lines[!grepl("^;", lines)]
    hdr   <- strsplit(trimws(lines[1]), "\\s+")[[1]]
    # RASP .eng header: name diam(mm) length(mm) delays prop_mass(kg) total_mass(kg) mfr
    # Per spec (thrustcurve.org/info/raspformat.html) mass fields are in KILOGRAMS.
    pm       <- as.numeric(hdr[5])          # propellant mass (kg)
    total_m  <- as.numeric(hdr[6])          # total motor mass — prop + casing (kg)
    casing_m <- max(total_m - pm, 0)        # casing: hardware that stays on board after burnout
    pairs <- lapply(lines[-1], function(l) as.numeric(strsplit(trimws(l), "\\s+")[[1]]))
    pairs <- Filter(function(p) length(p) >= 2 && !anyNA(p), pairs)
    list(
      thrust_curve   = data.frame(time = sapply(pairs, `[`, 1), thrust = sapply(pairs, `[`, 2)),
      prop_mass      = pm,
      casing_mass    = casing_m,
      motor_length_m = as.numeric(hdr[3]) / 1000,   # mm -> m
      motor_diam_m   = as.numeric(hdr[2]) / 1000    # mm -> m
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

# ── CSS ──────────────────────────────────────────────────────────────────────
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
"

# ── UI helper: paired numeric + unit selector ────────────────────────────────
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

# ── UI ───────────────────────────────────────────────────────────────────────
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
               column(5,
                      navset_card_pill(
                        nav_panel("Rocket",
                                  unit_input("dry_mass_val", "Dry mass",           90,  "g",  unit_choices_mass),
                                  unit_input("diameter",     "Body tube diameter", 24,  "mm", unit_choices_length),
                                  unit_input("body_length",  "Body tube length",   300, "mm", unit_choices_length),
                                  unit_input("cg_measured",  "CG from nose tip (unloaded)",   220, "mm", unit_choices_length)
                        ),
                        nav_panel("Nosecone",
                                  selectInput("nose_type","Nosecone type", choices=c("ogive","conical","parabolic")),
                                  unit_input("nose_length","Nosecone length", 70, "mm", unit_choices_length)
                        ),
                        nav_panel("Fins",
                                  numericInput("fin_count","Number of fins", value=3),
                                  unit_input("fin_root",  "Root chord",      50, "mm", unit_choices_length),
                                  unit_input("fin_tip",   "Tip chord",       25, "mm", unit_choices_length),
                                  unit_input("fin_span",  "Semi-span (from body wall to fin tip)", 30, "mm", unit_choices_length),
                                  unit_input("fin_sweep", "Sweep length (leading edge, axial projection)", 20, "mm", unit_choices_length)
                        ),
                        nav_panel("Chute",
                                  unit_input("parachute_diameter","Chute diameter", 305, "mm", unit_choices_length)
                        ),
                        nav_panel("Launch",
                                  unit_input("rail_length",    "Rail length",  0.9, "m",   unit_choices_length),
                                  unit_input("wind_speed_val", "Wind speed",   3,   "m/s", unit_choices_speed),
                                  sliderInput("wind_dir","Wind from (° CW from N)", min=0, max=360, value=270),
                                  numericInput("launch_angle","Launch angle from vertical (°)", value=0, min=0, max=30),
                                  sliderInput("launch_bearing","Launch bearing (° CW from N)", value=0, min=0, max=360)
                        ),
                        nav_panel("Engine",
                                  fileInput("motor_file", NULL, accept=".eng", buttonLabel="Upload .eng"),
                                  numericInput("parachute_delay","Ejection delay (s)", value=4),
                                  selectInput("engine_choice","Or choose an engine:",
                                              choices=c("Select engine..."="","A8","A10","B4","B6","C6","C11",
                                                        "D12","E12","E16","G40"),
                                              selected="B6", size=6, selectize=FALSE)
                        )
                      ),
                      br(),
                      uiOutput("stability_indicator")
               ),
               column(7,
                      plotOutput("thrust_curve_plot", height="400px"),
                      br(),
                      tags$img(src="fin_diagram.jpg", style="width:60%;border-radius:6px;border:1px solid #fc913a44;"),
                      p(style="color:#eae374aa;font-size:0.65rem;margin-top:6px;",
                        tags$a(href="https://www.apogeerockets.com/Peak-of-Flight/Newsletter615",
                               target="_blank", style="color:#D3D3D3;", "John K. Bennett; Apogee Peak of Flight Issue #615")
                      )
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
                          sliderInput("precision","Integration interval (s)", value=0.05, min=0.01, max=0.1),
                          p("0.05 s recommended"),
                          div(style="margin:10px 0 6px;", tags$label(class="unit-lbl","Display units")),
                          radioButtons("units", label=NULL,
                                       choices=c("Metric"="metric","Imperial (ft)"="imperial"),
                                       selected="metric", inline=FALSE),
                          br(),
                          actionButton("run","> Simulate", class="btn-primary", style="width:100%;")
                      ),
                      br(),
                      conditionalPanel("output.has_results",
                                       div(class="card", style="padding:16px;",
                                           h6(style="color:var(--c3);text-transform:uppercase;letter-spacing:1px;font-size:0.72rem;margin-bottom:10px;",
                                              "Flight summary"),
                                           verbatimTextOutput("summary")
                                       )
                      )
               ),
               column(9,
                      conditionalPanel("output.has_results",
                                       fluidRow(
                                         column(6, plotOutput("altitude_plot", height="260px")),
                                         column(6, plotOutput("velocity_plot", height="260px"))
                                       ),
                                       br(),
                                       plotlyOutput("track_3d", height="420px")
                      ),
                      conditionalPanel("!output.has_results",
                                       div(style=paste0(
                                         "display:flex;align-items:center;justify-content:center;",
                                         "height:500px;color:#fc913a44;font-size:0.85rem;",
                                         "text-transform:uppercase;letter-spacing:2px;"
                                       ), "▶  Press Simulate to run a flight")
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
                 p("0.05 s step recommended for speed/accuracy balance."),
                 sliderInput("chute_delay_std_dev",  "Ejection delay sd (s)",    value=1,  min=0.1, max=5),
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
                       "RRRocket 3D simulates model rocket flights, accounting for aerodynamic drag,
            motor thrust curves, wind weathercocking, parachute descent, and atmospheric
            density decreasing with altitude. A Monte Carlo framework quantifies uncertainty
            by varying launch conditions to produce a probabilistic landing footprint."),
                     br(),
                     h6(style="color:var(--c3);text-transform:uppercase;letter-spacing:1px;font-size:0.72rem;","Physics"),
                     p(style="color:var(--text);font-size:0.85rem;line-height:1.8;",
                       "Aerodynamics use the standard Barrowman equations (1967). Drag is split into
            nosecone pressure drag, skin-friction drag (body, nosecone, fins — all
            computed from wetted area with Cf = 0.005), and base drag (Hoerner formula).
            Drag is Mach-corrected with Prandtl-Glauert below M=0.8 and an empirical
            transonic ramp above. Atmosphere follows the International Standard Atmosphere
            troposphere model. Wind turbulence uses an Ornstein-Uhlenbeck process with
            a 1/7-power-law wind-shear profile."),
                     br(),
                     h6(style="color:var(--c3);text-transform:uppercase;letter-spacing:1px;font-size:0.72rem;","Monte Carlo"),
                     p(style="color:var(--text);font-size:0.85rem;line-height:1.8;",
                       "Real launches have uncertainty. The Monte Carlo tool runs hundreds of
            simulations with small variations in launch conditions, then plots every
            projected landing point on a satellite map. Draw a polygon around your safe
            landing zone and see what percentage of flights land inside it."),
                     br(),
                     h6(style="color:var(--c3);text-transform:uppercase;letter-spacing:1px;font-size:0.72rem;","Engine data"),
                     p(style="color:var(--text);font-size:0.85rem;line-height:1.8;",
                       "Motors use the standard RASP .eng format. Built-in curves cover common Estes
            A–G motors. For anything larger, download the .eng file from ",
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
             ), "<GitHub>")
    )
  )
)

# ── Plot theme ───────────────────────────────────────────────────────────────
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
    plot.margin      = margin(8,12,8,8))
}

# ── Server ───────────────────────────────────────────────────────────────────
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
  
  # ── SI stores (all values in SI units internally) ──────────────────────────
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
    parachute_diam = reactiveVal(0.305),
    rail_length    = reactiveVal(0.900),
    wind_speed     = reactiveVal(3.000)
  )
  
  # Generic observer factories — update SI store and convert display when unit changes
  make_length_observer <- function(id, store) {
    prev_unit <- reactiveVal(NULL)
    observeEvent(input[[paste0(id,"_unit")]], {
      pu <- prev_unit(); nu <- input[[paste0(id,"_unit")]]
      if (!is.null(pu) && isTruthy(input[[id]])) {
        si_val <- to_meters(input[[id]], pu); store(si_val)
        updateNumericInput(session, id, value=round(from_meters(si_val, nu), 4))
      }
      prev_unit(nu)
    }, ignoreInit=FALSE)
    observeEvent(input[[id]], {
      u <- input[[paste0(id,"_unit")]]
      if (isTruthy(u) && isTruthy(input[[id]])) store(to_meters(input[[id]], u))
    }, ignoreInit=TRUE)
  }
  
  make_mass_observer <- function(id, store) {
    prev_unit <- reactiveVal(NULL)
    observeEvent(input[[paste0(id,"_unit")]], {
      pu <- prev_unit(); nu <- input[[paste0(id,"_unit")]]
      if (!is.null(pu) && isTruthy(input[[id]])) {
        si_val <- to_kg(input[[id]], pu); store(si_val)
        updateNumericInput(session, id, value=round(from_kg(si_val, nu), 4))
      }
      prev_unit(nu)
    }, ignoreInit=FALSE)
    observeEvent(input[[id]], {
      u <- input[[paste0(id,"_unit")]]
      if (isTruthy(u) && isTruthy(input[[id]])) store(to_kg(input[[id]], u))
    }, ignoreInit=TRUE)
  }
  
  make_speed_observer <- function(id, store) {
    prev_unit <- reactiveVal(NULL)
    observeEvent(input[[paste0(id,"_unit")]], {
      pu <- prev_unit(); nu <- input[[paste0(id,"_unit")]]
      if (!is.null(pu) && isTruthy(input[[id]])) {
        si_val <- to_ms(input[[id]], pu); store(si_val)
        updateNumericInput(session, id, value=round(from_ms(si_val, nu), 4))
      }
      prev_unit(nu)
    }, ignoreInit=FALSE)
    observeEvent(input[[id]], {
      u <- input[[paste0(id,"_unit")]]
      if (isTruthy(u) && isTruthy(input[[id]])) store(to_ms(input[[id]], u))
    }, ignoreInit=TRUE)
  }
  
  make_mass_observer(  "dry_mass_val",       si$dry_mass)
  make_length_observer("diameter",           si$diameter)
  make_length_observer("body_length",        si$body_length)
  make_length_observer("cg_measured",        si$cg_measured)
  make_length_observer("nose_length",        si$nose_length)
  make_length_observer("fin_root",           si$fin_root)
  make_length_observer("fin_tip",            si$fin_tip)
  make_length_observer("fin_span",           si$fin_span)
  make_length_observer("fin_sweep",          si$fin_sweep)
  make_length_observer("parachute_diameter", si$parachute_diam)
  make_length_observer("rail_length",        si$rail_length)
  make_speed_observer( "wind_speed_val",     si$wind_speed)
  
  # ── Motor data ─────────────────────────────────────────────────────────────
  motor_data <- reactive({
    if (!is.null(input$motor_file)) {
      tryCatch(parse_thrust_input(input$motor_file, NULL), error=function(e) NULL)
    } else {
      ec <- input$engine_choice
      if (is.null(ec) || ec == "") return(NULL)
      tryCatch(parse_thrust_input(NULL, ec), error=function(e) NULL)
    }
  })
  
  # ── Aero (dry rocket geometry) ─────────────────────────────────────────────
  aero_reactive <- reactive({
    if (!all(sapply(list(input$nose_type, input$fin_count), isTruthy))) return(NULL)
    if (!all(sapply(list(si$nose_length(), si$body_length(), si$diameter(),
                         si$fin_root(), si$fin_tip(), si$fin_span(),
                         si$fin_sweep(), si$cg_measured()),
                    function(x) isTruthy(x) && x > 0))) return(NULL)
    tryCatch(
      compute_aero(input$nose_type,
                   si$nose_length(), si$body_length(), si$diameter(), input$fin_count,
                   si$fin_root(), si$fin_tip(), si$fin_span(), si$fin_sweep(),
                   si$cg_measured()),
      error=function(e) NULL)
  })
  
  # ── Aero with loaded CG (motor installed) ──────────────────────────────────
  # The motor contributes two distinct mass components:
  #   • casing_mass  — hardware that is on board for the entire flight
  #   • prop_mass    — propellant that burns away; starts at full, ends at 0
  # Both sit at cg_motor (geometric center of the motor).
  # "Loaded" stability uses the full motor mass (casing + propellant) — worst case for
  # CP margin because this is the heaviest/most-aft CG the rocket will ever have.
  # "Burnout" stability uses casing only — often the most critical for stability
  # because prop mass is gone but CG has moved forward (lighter motor end).
  aero_loaded <- reactive({
    aero <- aero_reactive(); if (is.null(aero)) return(NULL)
    td   <- motor_data()
    if (is.null(td)) {
      return(c(aero, list(stability_margin_loaded  = aero$stability_margin,
                          stability_margin_burnout = aero$stability_margin,
                          cg_loaded                = si$cg_measured())))
    }
    cg_motor     <- si$nose_length() + si$body_length() - td$motor_length_m / 2
    motor_mass   <- td$prop_mass + td$casing_mass       # total motor mass at ignition
    # Loaded CG: airframe + full motor
    total_loaded <- si$dry_mass() + motor_mass
    cg_loaded    <- (si$dry_mass() * si$cg_measured() + motor_mass * cg_motor) / total_loaded
    # Burnout CG: airframe + empty casing (propellant gone)
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
  
  # ── Stability indicator ────────────────────────────────────────────────────
  output$stability_indicator <- renderUI({
    al <- aero_loaded(); if (is.null(al)) return(NULL)
    sm  <- al$stability_margin_loaded
    fmt <- function(m) sprintf("%.1f mm  /  %.2f in", m*1000, m*39.3701)
    cls <- if (sm < 0.5) "stab-red" else if (sm < 1.0) "stab-yellow" else if (sm <= 3.0) "stab-green" else "stab-yellow"
    lbl <- if (sm < 0.5) "UNSTABLE" else if (sm < 1.0) "MARGINAL" else if (sm <= 3.0) "STABLE" else "OVERSTABLE"
    hint <- if (sm < 0.5) "Move CG forward or increase fin size." else
      if (sm > 3.0) "Risk of weathercocking in wind." else ""
    td <- motor_data()
    motor_lines <- if (!is.null(td)) {
      tagList(
        sprintf("CG loaded (ignition): %s", fmt(al$cg_loaded)), tags$br(),
        sprintf("Stability Margin when fully loaded:   %.2f cal", al$stability_margin_loaded), tags$br(),
        sprintf("Stability Margin at engine burnout:  %.2f cal", al$stability_margin_burnout), tags$br(),
        tags$span(style="color:var(--dim);font-size:0.72rem;",
                  sprintf("(casing %.0f g, prop %.0f g)",
                          td$casing_mass*1000, td$prop_mass*1000))
      )
    } else {
      tags$span(style="color:var(--dim);", "Load an engine to see loaded CG")
    }
    div(class=paste("stab-box", cls),
        tags$b(sprintf("%s — %.2f cal (loaded)", lbl, sm)), tags$br(),
        sprintf("CP: %s", fmt(al$CP)), tags$br(),
        motor_lines,
        if (nchar(hint) > 0) tagList(tags$br(), tags$span(hint)) else NULL
    )
  })
  
  # ── Thrust curve plot ──────────────────────────────────────────────────────
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
    t_ann <- if (use_metric()) sprintf("Total: %.2f Ns", ti) else sprintf("Total: %.2f lbf·s", ti*N_to_lbf)
    p_ann <- if (use_metric()) sprintf("Peak:  %.2f N",  mt) else sprintf("Peak:  %.2f lbf",   mt*N_to_lbf)
    ggplot(tc2, aes(time, thrust)) +
      geom_area(fill="#1a56db", alpha=0.08) +
      geom_line(color="#1a56db", linewidth=1) +
      geom_hline(yintercept=0, color="#e8eaf0") +
      annotate("text", x=max(tc2$time), y=max(tc2$thrust),
               label=t_ann, hjust=1, vjust=1.3, size=3.5, fontface="bold", color="#111928") +
      annotate("text", x=max(tc2$time), y=max(tc2$thrust)*0.87,
               label=p_ann, hjust=1, vjust=1.3, size=3.5, color="#6b7280") +
      scale_y_continuous(limits=c(0, max(tc2$thrust)*1.18)) +
      labs(x="time (s)", y=ylab, title="Engine thrust vs. time") +
      theme_plot()
  })
  
  # ── Simulate ───────────────────────────────────────────────────────────────
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
        aero$Cd_parasite, aero$CNa_total, aero$CP, input$precision,
        si$wind_speed(), input$wind_dir,
        si$cg_measured(), si$nose_length(), si$body_length(),
        parsed$motor_length_m, si$rail_length(),
        input$launch_bearing, input$launch_angle,
        landing_only=FALSE),
      error=function(e) { showNotification(paste("Sim error:", e$message), type="error"); NULL })
    if (is.null(sim)) return()
    res  <- list(sim=sim, aero=aero,
                 burn_time = max(parsed$thrust_curve$time),
                 label=paste0("Run ", length(run_history())+1),
                 motor=if (isTruthy(input$engine_choice) && input$engine_choice!="") input$engine_choice else "custom")
    hist <- run_history(); hist[[length(hist)+1]] <- res; run_history(hist)
    results_store(res)
  })
  
  # ── Monte Carlo ────────────────────────────────────────────────────────────
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
            parsed$prop_mass * rnorm(1, 1, 0.01*input$prop_mass_std_dev),
            si$dry_mass()    * rnorm(1, 1, 0.01*input$dry_mass_std_dev),
            parsed$casing_mass,   # casing mass is fixed — no uncertainty here
            si$diameter(), max(si$parachute_diam(), 0.05),
            max(0, rnorm(1, input$parachute_delay, input$chute_delay_std_dev)),
            aero$Cd_parasite * rnorm(1, 1, 0.01*input$cd_std_dev),
            aero$CNa_total, aero$CP, input$montecarlo_precision,
            max(0, rnorm(1, si$wind_speed(), 0.01*si$wind_speed()*input$wind_speed_std_dev)),
            rnorm(1, input$wind_dir, input$wind_dir_std_dev),
            si$cg_measured(), si$nose_length(), si$body_length(),
            parsed$motor_length_m, si$rail_length(),
            input$launch_bearing,
            max(0, rnorm(1, input$launch_angle, input$launch_angle_std_dev)),
            landing_only=TRUE),
          error=function(e) NULL)
        landings[[i]] <- if (!is.null(sim)) data.frame(x=sim$x, y=sim$y) else data.frame(x=0, y=0)
      }
    })
    mc_store(do.call(rbind, landings))
  })
  
  # ── Outputs ────────────────────────────────────────────────────────────────
  output$run_history_table <- renderUI({
    hist <- run_history(); if (length(hist) == 0) return(p("No runs yet."))
    sc <- if (use_metric()) 1 else m_to_ft
    u  <- if (use_metric()) "m" else "ft"
    rows <- lapply(rev(seq_along(hist)), function(i) {
      r <- hist[[i]]; s <- r$sim
      tags$tr(tags$td(r$label), tags$td(r$motor),
              tags$td(sprintf("%.0f %s", max(s$altitude)*sc, u)),
              tags$td(sprintf("%.1f s",  max(s$time))),
              tags$td(sprintf("%.2f cal", r$aero$stability_margin)))
    })
    tags$table(class="run-table",
               tags$thead(tags$tr(tags$th("Run"), tags$th("Motor"),
                                  tags$th("Apogee"), tags$th("Time"), tags$th("Stab."))),
               tags$tbody(rows))
  })
  
  output$summary <- renderPrint({
    res <- results_store(); req(!is.null(res))
    r <- res$sim; ae <- res$aero
    sc <- if (use_metric()) 1 else m_to_ft
    u  <- if (use_metric()) "m"   else "ft"
    us <- if (use_metric()) "m/s" else "ft/s"
    cat(sprintf("apogee             %d %s\n",       round(max(r$altitude)*sc), u))
    cat(sprintf("max velocity       %.1f %s\n",     max(r$velocity)*sc, us))
    rail_v <- attr(r, "rail_exit_ms")
    if (!is.null(rail_v) && !is.na(rail_v)) {
      flag <- if (rail_v*sc < (if (use_metric()) 15 else 49)) " *** LOW" else ""
      cat(sprintf("rail exit speed    %.1f %s%s\n", rail_v*sc, us, flag))
    }
    cat(sprintf("max Mach           %.3f\n",         max(r$mach)))
    cat(sprintf("time to apogee     %.2f s\n",       r$time[which.max(r$altitude)]))
    cat(sprintf("total flight time  %.2f s\n",       max(r$time)))
    # r$stability_margin is the live SM at each timestep; first value = at ignition (loaded),
    # last powered value ≈ burnout SM. min() gives the worst-case across the whole flight.
    cat(sprintf("SM at ignition     %.2f cal\n",     r$stability_margin[1]))
    cat(sprintf("SM at burnout      %.2f cal\n",     r$stability_margin[which.min(abs(r$time - max(r$time[r$time <= res$burn_time])))]))
    cat(sprintf("SM minimum         %.2f cal\n",     min(r$stability_margin)))
    cat(sprintf("Cd                 %.4f\n",         ae$Cd))
    cat(sprintf("  nose             %.4f\n",         ae$Cd_nose))
    cat(sprintf("  body             %.4f\n",         ae$Cd_body))
    cat(sprintf("  fins             %.4f\n",         ae$Cd_fins))
    cat(sprintf("  base             %.4f\n",         ae$Cd_base))
  })
  
  output$altitude_plot <- renderPlot({
    res <- results_store(); req(!is.null(res)); r <- res$sim
    alt  <- if (use_metric()) r$altitude else r$altitude * m_to_ft
    ylab <- if (use_metric()) "altitude (m)" else "altitude (ft)"
    ggplot(data.frame(t=r$time, alt=alt), aes(t, alt)) +
      geom_area(fill="#fc913a", alpha=0.15) +
      geom_line(color="#fc913a", linewidth=1) +
      geom_hline(yintercept=0, color="#e8eaf0") +
      labs(x="time (s)", y=ylab, title="Altitude") + theme_plot()
  })
  
  output$velocity_plot <- renderPlot({
    res <- results_store(); req(!is.null(res)); r <- res$sim
    vz   <- if (use_metric()) r$vz else r$vz * m_to_ft
    ylab <- if (use_metric()) "vertical velocity (m/s)" else "vertical velocity (ft/s)"
    ggplot(data.frame(t=r$time, vz=vz), aes(t, vz)) +
      geom_line(color="#f9d62e", linewidth=1) +
      geom_hline(yintercept=0, color="#fc913a44", linetype="dashed") +
      labs(x="time (s)", y=ylab, title="Vertical velocity") + theme_plot()
  })
  
  output$track_3d <- renderPlotly({
    res <- results_store(); req(!is.null(res)); r <- res$sim
    sc <- if (use_metric()) 1 else m_to_ft
    xl <- if (use_metric()) "East (m)"     else "East (ft)"
    yl <- if (use_metric()) "North (m)"    else "North (ft)"
    zl <- if (use_metric()) "Altitude (m)" else "Altitude (ft)"
    plot_ly(r, x=~x*sc, y=~y*sc, z=~altitude*sc, type="scatter3d", mode="lines",
            line=list(color=~altitude*sc,
                      colorscale=list(c(0,"#1a56db"), c(1,"#60a5fa")), width=3)) |>
      add_trace(x=tail(r$x,1)*sc, y=tail(r$y,1)*sc, z=0, type="scatter3d", mode="markers",
                marker=list(color="#e02424", size=5), name="landing") |>
      layout(paper_bgcolor="#ffffff", font=list(color="#111928"),
             scene=list(bgcolor="#f8f9fb",
                        xaxis=list(title=xl, gridcolor="#d0d5de"),
                        yaxis=list(title=yl, gridcolor="#d0d5de"),
                        zaxis=list(title=zl, gridcolor="#d0d5de")))
  })
  
  output$landing_pct <- renderPrint({
    lc <- mc_store(); req(!is.null(lc))
    sc <- if (use_metric()) 1 else m_to_ft
    u  <- if (use_metric()) "m" else "ft"
    drift <- sqrt(lc$x^2 + lc$y^2) * sc
    cat(sprintf("95th pct distance from pad: %.0f %s\n", quantile(drift, 0.95), u))
    cat(sprintf("Max distance from pad:      %.0f %s\n", max(drift), u))
    if (!is.null(landing_pct_val()))
      cat(sprintf("%.1f%% of rockets land in the polygon\n", landing_pct_val()))
  })
  
  # ── Map ────────────────────────────────────────────────────────────────────
  launch_point  <- reactiveVal(list(lat=38.89, lng=-77.03))
  drawn_polygon <- reactiveVal(NULL)
  
  observeEvent(input$map_click,         { launch_point(list(lat=input$map_click$lat, lng=input$map_click$lng)) })
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
    # Point-in-polygon (ray-casting):
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

#physics reference: Barrowman (1967), Hoerner (1965), ISA 1976
shinyApp(ui=ui, server=server)