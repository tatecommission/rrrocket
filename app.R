###NEW VERSION
library(shiny)
library(bslib)
library(plotly)
library(ggplot2)
library(leaflet)
library(leaflet.extras)

g0       <- 9.80665
R_air    <- 287.058
gamma_a  <- 1.4
m_to_ft  <- 3.28084
N_to_lbf <- 0.224809

unit_choices_length <- c("mm", "cm", "in", "ft", "m")
unit_choices_mass   <- c("g", "oz", "kg", "lb")
unit_choices_speed  <- c("m/s", "mph", "km/h", "knots")

to_meters <- function(val, unit) {
  if (!isTruthy(val)) return(0)
  switch(unit,
         "mm"    = val / 1000,
         "cm"    = val / 100,
         "in"    = val * 0.0254,
         "ft"    = val * 0.3048,
         "m"     = val,
         val / 1000
  )
}

to_kg <- function(val, unit) {
  if (!isTruthy(val)) return(0)
  switch(unit,
         "g"  = val / 1000,
         "oz" = val * 0.0283495,
         "kg" = val,
         "lb" = val * 0.453592,
         val / 1000
  )
}

to_ms <- function(val, unit) {
  if (!isTruthy(val)) return(0)
  switch(unit,
         "m/s"   = val,
         "mph"   = val * 0.44704,
         "km/h"  = val / 3.6,
         "knots" = val * 0.514444,
         val
  )
}

# Convert SI value back to display unit (for unit switching)
from_meters <- function(si_val, unit) {
  switch(unit,
         "mm" = si_val * 1000,
         "cm" = si_val * 100,
         "in" = si_val / 0.0254,
         "ft" = si_val / 0.3048,
         "m"  = si_val,
         si_val * 1000
  )
}

from_kg <- function(si_val, unit) {
  switch(unit,
         "g"  = si_val * 1000,
         "oz" = si_val / 0.0283495,
         "kg" = si_val,
         "lb" = si_val / 0.453592,
         si_val * 1000
  )
}

from_ms <- function(si_val, unit) {
  switch(unit,
         "m/s"   = si_val,
         "mph"   = si_val / 0.44704,
         "km/h"  = si_val * 3.6,
         "knots" = si_val / 0.514444,
         si_val
  )
}

isa_fast <- function(z) {
  z <- max(z, 0)
  if (z <= 11000) {
    T   <- 288.15 - 0.0065*z
    rho <- 1.225 * (T/288.15)^4.2561
  } else if (z <= 20000) {
    T   <- 216.65
    rho <- 0.36392 * exp(-0.0001577*(z - 11000))
  } else if (z <= 32000) {
    T   <- 216.65 + 0.001*(z - 20000)
    rho <- 0.08803 * (T/216.65)^(-34.1632/0.001)
  } else {
    T   <- 228.65 + 0.0028*(z - 32000)
    rho <- 0.01322 * (T/228.65)^(-34.1632/0.0028)
  }
  T   <- max(T, 150)
  rho <- max(rho, 1e-6)
  list(rho=rho, a=sqrt(gamma_a*R_air*T))
}

mach_cd_factor <- function(M) {
  ifelse(M<0.8, 1/sqrt(max(1-M^2,0.01)),
         ifelse(M<1.0, 1+1.4*(M-0.8)/0.2,
                ifelse(M<2.0, 2.4-0.8*(M-1.0), 1.6-0.2*(M-2.0))))
}

compute_aero <- function(nose_type, nose_length, body_length,
                         bt_diameter, fin_count, fin_root,
                         fin_tip, fin_span, fin_sweep, cg_measured) {
  bt_radius <- bt_diameter/2; Aref <- pi*bt_radius^2
  Xcp_nose <- switch(nose_type, conical=(2/3)*nose_length,
                     ogive=0.466*nose_length, parabolic=0.5*nose_length)
  ha <- atan(bt_radius/nose_length)
  Cd_nose <- switch(nose_type, conical=0.8*sin(ha)^2,
                    ogive=0.5*sin(ha)^2, parabolic=0.3*(bt_diameter/nose_length)^2)
  CNa_nb  <- 2.0
  CNa_bg  <- 2*0.9*(body_length*bt_diameter/Aref)
  Xcp_bg  <- nose_length + 0.45*body_length
  CNa_ne  <- CNa_nb + CNa_bg
  Xcp_ne  <- (CNa_nb*Xcp_nose + CNa_bg*Xcp_bg)/CNa_ne
  Xb      <- nose_length + body_length
  Afin    <- 0.5*(fin_root+fin_tip)*fin_span
  Kfb     <- 1+(bt_radius/(fin_span+bt_radius))
  CNa_fin <- Kfb*(4*fin_count*(fin_span/bt_diameter)^2)/
    (1+sqrt(1+(2*fin_span/(fin_root+fin_tip))^2))
  Xcp_fin <- Xb + (fin_sweep/3)*((fin_root+2*fin_tip)/(fin_root+fin_tip)) +
    (1/6)*(fin_root+fin_tip-(fin_root*fin_tip)/(fin_root+fin_tip))
  CNa_tot <- CNa_ne + CNa_fin
  Xcp     <- (CNa_ne*Xcp_ne + CNa_fin*Xcp_fin)/CNa_tot
  Cd_fins <- fin_count*2*Afin*0.008/Aref
  Cd_body <- 0.01
  Cd_base <- 0.029/sqrt(max(Cd_nose+Cd_body+Cd_fins,1e-6))
  list(Cd=Cd_nose+Cd_body+Cd_fins+Cd_base, CP=Xcp,
       stability_margin=(Xcp-cg_measured)/bt_diameter,
       CNa_total=CNa_tot,
       Cd_nose=Cd_nose, Cd_fins=Cd_fins, Cd_body=Cd_body, Cd_base=Cd_base)
}

weathercock_gain <- function(sm) max(min(sm/3.0,1.0)*0.6,0)

flight_simulation_3d <- function(thrust_curve, prop_mass, dry_mass,
                                 bt_diameter, chute_diameter, chute_delay,
                                 Cd, CNa_total, Xcp, precision,
                                 wind_speed_ref, wind_dir_deg,
                                 cg_dry_m, nose_length, body_length,
                                 motor_length_m, rail_length_m,
                                 launch_bearing_deg, launch_angle_deg,
                                 landing_only = FALSE) {

  if (is.null(thrust_curve) || nrow(thrust_curve) == 0) return(NULL)
  thrust    <- approxfun(thrust_curve$time, thrust_curve$thrust, yleft=0, yright=0)
  burn_time <- max(thrust_curve$time)
  m0        <- dry_mass + prop_mass
  bt_area   <- pi*(bt_diameter/2)^2
  chute_area <- pi*(chute_diameter/2)^2
  total_imp <- integrate(thrust, min(thrust_curve$time), max(thrust_curve$time))$value
  v_exhaust <- total_imp / prop_mass
  if (!is.finite(v_exhaust) || v_exhaust <= 0) {
    showNotification("check .eng file", type="error")
    return(NULL)
  }
  dt        <- max(precision, 0.005)
  ar <- launch_angle_deg*pi/180; br <- launch_bearing_deg*pi/180
  lx <- sin(ar)*sin(br); ly <- sin(ar)*cos(br); lz <- cos(ar)
  cg_motor <- nose_length + body_length - motor_length_m/2
  x<-0; y<-0; z<-0; vx<-0; vy<-0; vz<-0
  m <- m0; m_prop <- prop_mass; t <- 0; t_apogee <- NA; chute_open <- FALSE
  k_sh <- 2.0
  lam  <- max(wind_speed_ref, 0.01) / gamma(1 + 1/k_sh)
  weibull_mean <- lam * gamma(1 + 1/k_sh)
  wi   <- c(-wind_speed_ref*sin(wind_dir_deg*pi/180), -wind_speed_ref*cos(wind_dir_deg*pi/180))
  wx_s <- wi[1]; wy_s <- wi[2]
  on_rail  <- TRUE
  rail_ht  <- rail_length_m * cos(ar)
  max_steps <- ceiling(1200/dt)
  
  if (landing_only) {
    repeat {
      dt <- if (t <= burn_time + 0.5) min(precision, 0.01) else max(precision, 0.005)
      atm <- isa_fast(z); rho <- atm$rho; a_snd <- atm$a
      spd_ref <- wind_speed_ref * (max(z,10)/10)^0.14
      b_rad   <- wind_dir_deg*pi/180
      mu_x <- -spd_ref*sin(b_rad); mu_y <- -spd_ref*cos(b_rad)
      noise_amp <- rweibull(1, shape=k_sh, scale=lam) - weibull_mean
      wx_s <- wx_s + 0.5*(mu_x-wx_s)*dt + noise_amp*(mu_x/max(spd_ref,1e-6))*sqrt(dt)
      wy_s <- wy_s + 0.5*(mu_y-wy_s)*dt + noise_amp*(mu_y/max(spd_ref,1e-6))*sqrt(dt)
      vrx <- vx-wx_s; vry <- vy-wy_s; vrz <- vz
      vrm <- max(sqrt(vrx^2+vry^2+vrz^2), 1e-6); M <- vrm/a_snd
      Cd_eff   <- if(chute_open) 1.5 else Cd*mach_cd_factor(M)
      area_eff <- if(chute_open) chute_area else bt_area
      q   <- 0.5*rho*vrm^2; Fd <- Cd_eff*area_eff*q
      Fdx <- -Fd*(vrx/vrm); Fdy <- -Fd*(vry/vrm); Fdz <- -Fd*(vrz/vrm)
      Ft  <- thrust(t); Fgz <- -m*g0
      cg_live <- (dry_mass*cg_dry_m + m_prop*cg_motor)/m
      sm_live <- (Xcp-cg_live)/bt_diameter
      wc_gain <- weathercock_gain(sm_live)
      Fwx <- 0; Fwy <- 0
      if (on_rail && z >= rail_ht) on_rail <- FALSE
      if (!on_rail && is.na(t_apogee) && vz > 0.5) {
        lx2 <- wx_s-vx; ly2 <- wy_s-vy; lm2 <- max(sqrt(lx2^2+ly2^2), 1e-6)
        cp  <- wc_gain*q*bt_area*CNa_total
        Fwx <- cp*(lx2/lm2)*min(lm2/max(vrm,1), 1)
        Fwy <- cp*(ly2/lm2)*min(lm2/max(vrm,1), 1)
      }
      ax <- (Fdx+Fwx+Ft*lx)/m
      ay <- (Fdy+Fwy+Ft*ly)/m
      az <- (Ft*lz+Fdz+Fgz)/m
      if (t <= burn_time) { m_prop <- max(m_prop-(Ft/v_exhaust)*dt, 0); m <- dry_mass+m_prop }
      if (!is.na(t_apogee) && !chute_open && (t-t_apogee) >= chute_delay) chute_open <- TRUE
      vx <- vx+ax*dt; vy <- vy+ay*dt; vz <- vz+az*dt
      x  <- x+vx*dt;  y  <- y+vy*dt;  z  <- z+vz*dt; t <- t+dt
      if (is.na(t_apogee) && z > 1 && vz < 0) t_apogee <- t
      if (!is.na(t_apogee) && z <= 0) {
        frac <- if (abs(vz) > 1e-6) (z / (vz*dt)) else 0
        x <- x - vx*dt*frac
        y <- y - vy*dt*frac
        break
      }
      if (t > 1200) break
    }
    return(data.frame(x=x, y=y))
  }
  
  out <- data.frame(time=numeric(max_steps), x=numeric(max_steps), y=numeric(max_steps),
                    altitude=numeric(max_steps), velocity=numeric(max_steps),
                    vx=numeric(max_steps), vy=numeric(max_steps), vz=numeric(max_steps),
                    mach=numeric(max_steps), stability_margin=numeric(max_steps))
  i <- 1
  repeat {
    dt <- if (t <= burn_time + 0.5) min(precision, 0.01) else max(precision, 0.005)
    atm <- isa_fast(z); rho <- atm$rho; a_snd <- atm$a
    spd_ref <- wind_speed_ref*(max(z,10)/10)^0.14
    b_rad   <- wind_dir_deg*pi/180
    mu_x <- -spd_ref*sin(b_rad); mu_y <- -spd_ref*cos(b_rad)
    noise_amp <- rweibull(1, shape=k_sh, scale=lam) - weibull_mean
    wx_s <- wx_s + 0.5*(mu_x-wx_s)*dt + noise_amp*(mu_x/max(spd_ref,1e-6))*sqrt(dt)
    wy_s <- wy_s + 0.5*(mu_y-wy_s)*dt + noise_amp*(mu_y/max(spd_ref,1e-6))*sqrt(dt)
    vrx <- vx-wx_s; vry <- vy-wy_s; vrz <- vz
    vrm <- max(sqrt(vrx^2+vry^2+vrz^2), 1e-6); M <- vrm/a_snd
    Cd_eff   <- if(chute_open) 1.5 else Cd*mach_cd_factor(M)
    area_eff <- if(chute_open) chute_area else bt_area
    q   <- 0.5*rho*vrm^2; Fd <- Cd_eff*area_eff*q
    Fdx <- -Fd*(vrx/vrm); Fdy <- -Fd*(vry/vrm); Fdz <- -Fd*(vrz/vrm)
    Ft  <- thrust(t); Fgz <- -m*g0
    cg_live <- (dry_mass*cg_dry_m + m_prop*cg_motor)/m
    sm_live <- (Xcp-cg_live)/bt_diameter
    wc_gain <- weathercock_gain(sm_live)
    Fwx <- 0; Fwy <- 0
    if (on_rail && z >= rail_ht) on_rail <- FALSE
    if (!on_rail && is.na(t_apogee) && vz > 0.5) {
      lx2 <- wx_s-vx; ly2 <- wy_s-vy; lm2 <- max(sqrt(lx2^2+ly2^2), 1e-6)
      cp  <- wc_gain*q*bt_area*CNa_total
      Fwx <- cp*(lx2/lm2)*min(lm2/max(vrm,1), 1)
      Fwy <- cp*(ly2/lm2)*min(lm2/max(vrm,1), 1)
    }
    ax <- (Fdx+Fwx+Ft*lx)/m; ay <- (Fdy+Fwy+Ft*ly)/m; az <- (Ft*lz+Fdz+Fgz)/m
    out$time[i]<-t; out$x[i]<-x; out$y[i]<-y; out$altitude[i]<-z
    out$velocity[i]<-sqrt(vx^2+vy^2+vz^2)
    out$vx[i]<-vx; out$vy[i]<-vy; out$vz[i]<-vz
    out$mach[i]<-M; out$stability_margin[i]<-sm_live
    i <- i+1
    if (t <= burn_time) { m_prop <- max(m_prop-(Ft/v_exhaust)*dt, 0); m <- dry_mass+m_prop }
    if (!is.na(t_apogee) && !chute_open && (t-t_apogee) >= chute_delay) chute_open <- TRUE
    vx <- vx+ax*dt; vy <- vy+ay*dt; vz <- vz+az*dt
    x  <- x+vx*dt;  y  <- y+vy*dt;  z  <- z+vz*dt; t <- t+dt
    if (is.na(t_apogee) && z > 1 && vz < 0) t_apogee <- t
    if (!is.na(t_apogee) && z <= 0) break
    if (t > 1200 || i >= max_steps) break
  }
  result <- out[1:(i-1),]
  rail_row <- which(result$altitude >= rail_ht)
  attr(result,"rail_exit_ms") <- if(length(rail_row)>0) result$velocity[rail_row[1]] else NA_real_
  result
}

parse_thrust_input <- function(motor_file, engine_choice) {
  read_eng <- function(lines) {
    lines <- lines[!grepl("^;",lines)]
    hdr   <- strsplit(trimws(lines[1]),"\\s+")[[1]]
    pm <- as.numeric(hdr[5]) / 1000
    pairs <- lapply(lines[-1],function(l) as.numeric(strsplit(trimws(l),"\\s+")[[1]]))
    pairs <- Filter(function(p) length(p) >= 2 && !anyNA(p), pairs)
    list(thrust_curve=data.frame(time=sapply(pairs,`[`,1),thrust=sapply(pairs,`[`,2)),
         prop_mass=pm,
         motor_length_m=as.numeric(hdr[3])/1000,
         motor_diam_m=as.numeric(hdr[2])/1000)
  }
  if (!is.null(motor_file)) {
    read_eng(readLines(motor_file$datapath,warn=FALSE))
  } else {
    if (is.null(engine_choice) || engine_choice == "") return(NULL)
    path <- file.path("data",paste0(engine_choice,".eng"))
    if (!file.exists(path)) return(NULL)
    read_eng(readLines(path,warn=FALSE))
  }
}

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
      #ff4e50 0%,
      #ff4e50 16%,
      #fc913a 16%,
      #fc913a 32%,
      #f9d62e 32%,
      #f9d62e 48%,
      #eae374 48%,
      #eae374 64%,
      #e2f4c7 64%,
      #e2f4c7 80%,
      #b8dba0 80%,
      #b8dba0 100%) fixed!important;
}

body {
  background:transparent!important;
  color:var(--text)!important;
  min-height:100vh;
}

/* Faint inner-shadow on each band to give 3-D depth */
body::before {
  content:'';
  position:fixed;
  inset:0;
  background:
    repeating-linear-gradient(
      180deg,
      transparent 0px,
      transparent calc(16vh - 4px),
      rgba(0,0,0,0.22) calc(16vh - 4px),
      rgba(255,255,255,0.06) calc(16vh),
      transparent calc(16vh + 1px)
    );
  pointer-events:none;
  z-index:0;
}

/* overlay so content stays readable */
body::after {
  content:'';
  position:fixed;
  inset:0;
  background:rgba(10,6,2,0.55);
  pointer-events:none;
  z-index:1;
}

.tab-content,
.navbar,
.shiny-plot-output,
.card,
.well,
.leaflet-container,
pre,
.shiny-verbatim-output,
.stab-box,
.run-table,
.irs--shiny,
.form-control,
input,select,
.nav-pills,
.sidebarPanel,
.mainPanel,
.container-fluid,
#shiny-tab-Setup,
#shiny-tab-Simulation,
#shiny-tab-Results,
#shiny-tab-3D.Track,
#shiny-tab-Monte.Carlo {
  position:relative;
  z-index:2;
}

/* ── NAVBAR ── */
.navbar {
  background:linear-gradient(180deg,rgba(26,16,8,0.97) 0%,rgba(20,10,4,0.99) 100%)!important;
  border-bottom:2px solid var(--c2)!important;
  box-shadow:0 4px 24px rgba(255,78,80,0.25),0 2px 8px rgba(0,0,0,0.7)!important;
  padding:4px 32px!important;
  z-index:100!important;
  position:relative!important;
  min-height:58px!important;
}
.navbar-brand {
  display:flex!important;
  align-items:center!important;
  padding:8px 0!important;
  gap:0!important;
}
.nav-link {
  font-size:0.73rem!important;
  color:#f9d62eaa!important;
  padding:16px 20px!important;
  text-transform:uppercase!important;
  letter-spacing:1.4px!important;
  font-weight:600!important;
  border-bottom:3px solid transparent!important;
  transition:all .2s!important;
}
.nav-link:hover  { color:var(--c3)!important; }
.nav-link.active {
  color:var(--c1)!important;
  border-bottom-color:var(--c1)!important;
  background:transparent!important;
  text-shadow:0 0 12px rgba(255,78,80,0.5)!important;
}

/* ── SKEUOMORPHIC BUTTONS ── */
.btn {
  font-size:0.73rem!important;
  text-transform:uppercase!important;
  letter-spacing:1.1px!important;
  font-weight:700!important;
  border-radius:7px!important;
  padding:10px 26px!important;
  transition:all .08s ease!important;
  cursor:pointer!important;
}

/* Primary — sunset red */
.btn-primary {
  background:linear-gradient(180deg,#ff6b5b 0%,#ff4e50 50%,#d93535 100%)!important;
  border:1px solid #b02020!important;
  color:#fff!important;
  box-shadow:
    0 1px 0 rgba(255,255,255,0.3) inset,
    0 -2px 0 rgba(0,0,0,0.3) inset,
    0 4px 10px rgba(255,78,80,0.5),
    0 2px 4px rgba(0,0,0,0.5)!important;
  text-shadow:0 -1px 0 rgba(0,0,0,0.4)!important;
}
.btn-primary:hover {
  background:linear-gradient(180deg,#ff8070 0%,#ff6b5b 50%,#ff4e50 100%)!important;
  box-shadow:
    0 1px 0 rgba(255,255,255,0.35) inset,
    0 -2px 0 rgba(0,0,0,0.25) inset,
    0 6px 16px rgba(255,78,80,0.6),
    0 2px 6px rgba(0,0,0,0.4)!important;
}
.btn-primary:active {
  background:linear-gradient(180deg,#d93535 0%,#b02020 100%)!important;
  box-shadow:0 2px 4px rgba(0,0,0,0.6) inset!important;
  transform:translateY(1px)!important;
}

/* Warning — orange */
.btn-warning {
  background:linear-gradient(180deg,#ffaa55 0%,#fc913a 50%,#d97020 100%)!important;
  border:1px solid #b05810!important;
  color:#fff!important;
  box-shadow:
    0 1px 0 rgba(255,255,255,0.3) inset,
    0 -2px 0 rgba(0,0,0,0.3) inset,
    0 4px 10px rgba(252,145,58,0.5),
    0 2px 4px rgba(0,0,0,0.5)!important;
  text-shadow:0 -1px 0 rgba(0,0,0,0.35)!important;
}
.btn-warning:hover {
  background:linear-gradient(180deg,#ffc070 0%,#ffaa55 50%,#fc913a 100%)!important;
  box-shadow:
    0 1px 0 rgba(255,255,255,0.35) inset,
    0 -2px 0 rgba(0,0,0,0.25) inset,
    0 6px 16px rgba(252,145,58,0.55),
    0 2px 6px rgba(0,0,0,0.4)!important;
}
.btn-warning:active {
  background:linear-gradient(180deg,#d97020 0%,#b05810 100%)!important;
  box-shadow:0 2px 4px rgba(0,0,0,0.6) inset!important;
  transform:translateY(1px)!important;
}

/* Default */
.btn-default,.btn-secondary {
  background:linear-gradient(180deg,#3a2510 0%,#251508 100%)!important;
  border:1px solid #fc913a66!important;
  color:#f9d62e!important;
  box-shadow:
    0 1px 0 rgba(255,255,255,0.08) inset,
    0 -1px 0 rgba(0,0,0,0.4) inset,
    0 2px 6px rgba(0,0,0,0.5)!important;
}
.btn-default:hover,.btn-secondary:hover {
  background:linear-gradient(180deg,#4a3018 0%,#3a2510 100%)!important;
  color:var(--c3)!important;
}
.btn-default:active,.btn-secondary:active {
  box-shadow:0 2px 4px rgba(0,0,0,0.6) inset!important;
  transform:translateY(1px)!important;
}
.btn-file {
  background:linear-gradient(180deg,#3a2510,#251508)!important;
  border:1px solid #fc913a55!important;
  color:#f9d62e!important;
  box-shadow:0 2px 6px rgba(0,0,0,0.5)!important;
}

/* ── INPUTS ── */
label,.form-label {
  font-size:0.68rem!important;
  color:var(--c4)!important;
  text-transform:uppercase!important;
  letter-spacing:.08em!important;
  font-weight:600!important;
}
.form-control,.form-select,input[type=number],select {
  background:rgba(20,10,4,0.85)!important;
  border:1px solid #fc913a55!important;
  color:var(--text)!important;
  border-radius:5px!important;
  padding:6px 10px!important;
  font-size:0.84rem!important;
  box-shadow:0 2px 5px rgba(0,0,0,0.5) inset,0 1px 0 rgba(255,255,255,0.04)!important;
  transition:border-color .15s,box-shadow .15s!important;
}
.form-control:focus,input:focus,select:focus {
  border-color:var(--c2)!important;
  box-shadow:0 2px 5px rgba(0,0,0,0.5) inset,0 0 0 3px rgba(252,145,58,0.25)!important;
  outline:none!important;
}
select option { background:#1a0f05; color:var(--text); }

/* ── SLIDERS ── */
.irs--shiny .irs-bar         { background:linear-gradient(90deg,var(--c1),var(--c2))!important; border-color:var(--c1)!important; }
.irs--shiny .irs-handle      { border-color:var(--c2)!important; background:#fff8f0!important; box-shadow:0 2px 6px rgba(0,0,0,0.5)!important; }
.irs--shiny .irs-single      { background:var(--c2)!important; color:#1a1008!important; font-size:0.65rem!important; font-weight:700!important; }
.irs--shiny .irs-line        { background:#2a1508!important; box-shadow:0 1px 4px rgba(0,0,0,0.5) inset!important; }
.irs--shiny .irs-grid-text,
.irs--shiny .irs-min,
.irs--shiny .irs-max         { color:var(--c4)!important; background:transparent!important; }
.irs--shiny .irs-grid-pol    { background:var(--c4)!important; }

/* ── CARDS / PANELS ── */
.card {
  background:linear-gradient(160deg,rgba(30,15,5,0.92) 0%,rgba(18,9,3,0.95) 100%)!important;
  border:1px solid #fc913a44!important;
  border-radius:8px!important;
  box-shadow:0 6px 20px rgba(0,0,0,0.6),0 1px 0 rgba(255,200,100,0.06) inset!important;
}
.card-body { padding:16px!important; }
.well {
  background:rgba(18,9,3,0.88)!important;
  border:1px solid #fc913a44!important;
  border-radius:8px!important;
  box-shadow:0 3px 10px rgba(0,0,0,0.5) inset!important;
}

/* ── NAV PILLS ── */
.nav-pills .nav-link {
  font-size:0.67rem!important;
  color:var(--c4)!important;
  border-radius:5px!important;
  padding:5px 13px!important;
  text-transform:uppercase!important;
  letter-spacing:.08em!important;
  background:linear-gradient(180deg,#3a2010 0%,#1e0e04 100%)!important;
  border:1px solid #fc913a44!important;
  margin:2px!important;
  font-weight:600!important;
  box-shadow:0 2px 5px rgba(0,0,0,0.4),0 1px 0 rgba(255,200,80,0.07) inset!important;
  transition:all .12s!important;
}
.nav-pills .nav-link:hover {
  background:linear-gradient(180deg,#4a2a14 0%,#2e1608 100%)!important;
  color:var(--c3)!important;
}
.nav-pills .nav-link.active {
  background:linear-gradient(180deg,#ff6b5b 0%,#ff4e50 50%,#d93535 100%)!important;
  border-color:#b02020!important;
  color:#fff!important;
  box-shadow:0 3px 10px rgba(255,78,80,0.45),0 1px 0 rgba(255,255,255,0.2) inset!important;
}

/* ── TEXT ── */
p { font-size:0.7rem; color:var(--c4); margin:2px 0 8px; }
.tab-content { padding:20px 24px; }

/* ── VERBATIM OUTPUT ── */
pre,.shiny-verbatim-output {
  background:rgba(8,4,0,0.92)!important;
  color:#f9d62e!important;
  border:1px solid #fc913a44!important;
  border-radius:6px!important;
  font-family:'Courier New',monospace!important;
  font-size:0.78rem!important;
  padding:14px!important;
  line-height:1.9!important;
  box-shadow:0 3px 10px rgba(0,0,0,0.6) inset!important;
}

.shiny-plot-output {
  border:1px solid #fc913a44;
  border-radius:6px;
  box-shadow:0 6px 16px rgba(0,0,0,0.5);
}
.leaflet-container {
  border-radius:6px;
  border:1px solid #fc913a55;
  box-shadow:0 6px 16px rgba(0,0,0,0.5);
}

/* ── HOME ── */
.home-wrap {
  display:flex; flex-direction:column;
  justify-content:center; align-items:center;
  min-height:80vh; text-align:center;
  position:relative;
}

/* sunset graphic behind title */
.home-graphic {
  position:absolute;
  width:520px; height:520px;
  border-radius:50%;
  background:radial-gradient(ellipse at 50% 60%,
    rgba(249,214,46,0.18) 0%,
    rgba(252,145,58,0.14) 35%,
    rgba(255,78,80,0.10) 60%,
    transparent 75%);
  box-shadow:
    0 0 80px 20px rgba(249,214,46,0.12),
    0 0 160px 60px rgba(252,145,58,0.08),
    0 0 260px 100px rgba(255,78,80,0.05);
  top:50%; left:50%;
  transform:translate(-50%,-50%);
  pointer-events:none;
  z-index:0;
}

/* horizon line */
.home-horizon {
  position:absolute;
  bottom:28%; left:0; right:0;
  height:1px;
  background:linear-gradient(90deg,
    transparent 0%, rgba(249,214,46,0.35) 20%,
    rgba(252,145,58,0.6) 50%,
    rgba(249,214,46,0.35) 80%, transparent 100%);
  box-shadow:0 0 12px 2px rgba(252,145,58,0.3);
  pointer-events:none;
}

/* silhouette mountains */
.home-mountains {
  position:absolute;
  bottom:28%; left:0; right:0;
  height:80px;
  background:
    polygon(0% 100%,8% 40%,16% 65%,26% 20%,36% 55%,
            46% 30%,56% 60%,66% 15%,76% 50%,86% 35%,100% 55%,100% 100%);
  clip-path:polygon(
    0% 100%,  6% 55%,  12% 75%,
    20% 30%,  28% 60%,  36% 20%,
    44% 50%,  52% 35%,  60% 55%,
    68% 10%,  76% 45%,  84% 30%,
    92% 50%,  100% 40%, 100% 100%
  );
  background:#1a0a02;
  opacity:0.85;
  pointer-events:none;
  z-index:1;
}

.home-content { position:relative; z-index:2; }
.title-wrap { display:flex; align-items:center; justify-content:center; }
.home-title {
  font-size:5rem;
  font-weight:800;
  letter-spacing:3px;
  margin:0;
  background:linear-gradient(135deg,#ff4e50 0%,#fc913a 35%,#f9d62e 65%,#e2f4c7 100%);
  -webkit-background-clip:text;
  -webkit-text-fill-color:transparent;
  text-rendering:geometricPrecision;
  filter:drop-shadow(0 0 20px rgba(252,145,58,0.4));
}
.cursor {
  display:inline-block; width:4px; height:5rem;
  background:linear-gradient(180deg,#f9d62e,#fc913a);
  margin-left:8px; vertical-align:middle;
  animation:blink .7s step-end infinite;
  box-shadow:0 0 10px rgba(249,214,46,0.6);
  border-radius:2px;
}
@keyframes blink { 50% { opacity:0; } }
.home-sub {
  font-size:0.85rem; letter-spacing:4px;
  margin-top:12px; text-transform:uppercase; font-weight:500;
  color:var(--c4);
  text-shadow:0 0 20px rgba(234,227,116,0.4);
}

/* ── STABILITY BOX ── */
.stab-box {
  border-radius:6px; padding:12px 16px; margin-top:12px;
  font-size:0.78rem; line-height:2; font-weight:500;
  box-shadow:0 4px 14px rgba(0,0,0,0.5);
}
.stab-green  { color:#86efac; border:1px solid #22c55e; background:linear-gradient(160deg,rgba(13,40,24,0.9),rgba(19,34,24,0.92)); }
.stab-yellow { color:#fde68a; border:1px solid #f59e0b; background:linear-gradient(160deg,rgba(43,31,8,0.9),rgba(43,33,17,0.92)); }
.stab-red    { color:#fca5a5; border:1px solid #ef4444; background:linear-gradient(160deg,rgba(42,14,14,0.9),rgba(42,21,21,0.92)); }

/* ── RUN TABLE ── */
.run-table { width:100%; border-collapse:collapse; font-size:0.73rem; margin-top:10px; }
.run-table th {
  color:var(--c4); text-transform:uppercase; letter-spacing:.07em;
  padding:6px 10px; border-bottom:2px solid #fc913a55; text-align:left; font-weight:700;
}
.run-table td {
  color:var(--text); padding:6px 10px;
  border-bottom:1px solid rgba(252,145,58,0.15); font-weight:500;
}
.run-table tr:hover td { background:rgba(252,145,58,0.07); }

/* ── UNIT INPUT ── */
.unit-row-wrap { margin-bottom:10px; }
.unit-row { display:flex; gap:6px; align-items:center; }
.num-wrap { flex:2; }
.sel-wrap { flex:1; }
.unit-lbl {
  font-size:0.67rem; color:var(--c4); text-transform:uppercase;
  letter-spacing:.08em; font-weight:600; display:block; margin-bottom:3px;
}

/* stop run summary table from bleeding onto right */
.run-table { table-layout:fixed; word-break:break-word; }
.run-table td, .run-table th { overflow:hidden; text-overflow:ellipsis; white-space:nowrap; }

.navbar-nav.navbar-right-custom {
  margin-left: auto !important;
}
"

unit_input <- function(input_id, label, default_val, default_unit, choices) {
  div(class="unit-row-wrap",
      tags$label(class="unit-lbl", label),
      div(class="unit-row",
          div(class="num-wrap",
              numericInput(input_id, label=NULL, value=default_val, width="100%")),
          div(class="sel-wrap",
              selectInput(paste0(input_id,"_unit"), label=NULL,
                          choices=choices, selected=default_unit, width="100%"))
      )
  )
}

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
          var i = 0;
          el.textContent = '';
          function type() {
            if (i < text.length) {
              el.textContent += text[i++];
              setTimeout(type, 100);
            }
          }
          setTimeout(type, 300);
        });
      ")),
             div(class="home-wrap",
                 div(class="title-wrap",
                     h1(class="home-title", span(id="typed-title")),
                     span(class="cursor")
                 ),
                 p(class="home-sub","Model rocket flight simulator")
             )
    ),
    
    tabPanel("Setup",
             fluidRow(
               column(5,
                      navset_card_pill(
                        nav_panel("Rocket",
                                  unit_input("dry_mass_val", "Dry mass",            90,  "g",   unit_choices_mass),
                                  unit_input("diameter",     "Body tube diameter",  24,  "mm",  unit_choices_length),
                                  unit_input("body_length",  "Body tube length",    300, "mm",  unit_choices_length),
                                  unit_input("cg_measured",  "CG from nose tip",    220, "mm",  unit_choices_length)
                        ),
                        nav_panel("Nosecone",
                                  selectInput("nose_type","Nosecone type",choices=c("ogive","conical","parabolic")),
                                  unit_input("nose_length","Nosecone length",70,"mm",unit_choices_length)
                        ),
                        nav_panel("Fins",
                                  numericInput("fin_count","Number of fins",value=3),
                                  unit_input("fin_root",  "Root chord",   50, "mm", unit_choices_length),
                                  unit_input("fin_tip",   "Tip chord",    25, "mm", unit_choices_length),
                                  unit_input("fin_span",  "Semi-span",    30, "mm", unit_choices_length),
                                  unit_input("fin_sweep", "Sweep length", 20, "mm", unit_choices_length)
                        ),
                        nav_panel("Chute",
                                  unit_input("parachute_diameter","Chute diameter",305,"mm",unit_choices_length)
                        ),
                        nav_panel("Launch",
                                  unit_input("rail_length",   "Rail length", 0.9, "m",   unit_choices_length),
                                  unit_input("wind_speed_val","Wind speed",  3,   "m/s", unit_choices_speed),
                                  sliderInput("wind_dir","Wind from (° CW from N)",min=0,max=360,value=270),
                                  numericInput("launch_angle","Launch angle from vertical (°)",value=0,min=0,max=30),
                                  sliderInput("launch_bearing","Launch bearing (° CW from N)",value=0,min=0,max=360)
                        ),
                        nav_panel("Engine",
                                  fileInput("motor_file",NULL,accept=".eng",buttonLabel="Upload .eng"),
                                  numericInput("parachute_delay","Ejection delay (s)",value=4),
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
                      plotOutput("thrust_curve_plot", height="400px")
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
                          sliderInput("precision","Integration interval (s)",value=0.01,min=0.001,max=0.1),
                          p("0.01 s recommended"),
                          div(style="margin:10px 0 6px;",
                              tags$label(class="unit-lbl","Display units")),
                          radioButtons("units",label=NULL,
                                       choices=c("Metric"="metric","Imperial (ft)"="imperial"),
                                       selected="metric",inline=FALSE),
                          br(),
                          actionButton("run","> Simulate",class="btn-primary",style="width:100%;")
                      ),
                      br(),
                      conditionalPanel("output.has_results",
                                       div(class="card",style="padding:16px;",
                                           h6(style="color:var(--c3);text-transform:uppercase;letter-spacing:1px;font-size:0.72rem;margin-bottom:10px;",
                                              "Flight summary"),
                                           verbatimTextOutput("summary")
                                       )
                      )
               ),
               column(9,
                      conditionalPanel("output.has_results",
                                       fluidRow(
                                         column(6, plotOutput("altitude_plot",height="260px")),
                                         column(6, plotOutput("velocity_plot",height="260px"))
                                       ),
                                       br(),
                                       plotlyOutput("track_3d",height="420px")
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
                              div(class="card",style="padding:16px;overflow-x:auto;",
                                  h6(style="color:var(--c3);text-transform:uppercase;letter-spacing:1px;font-size:0.72rem;margin-bottom:10px;",
                                     "Run history"),
                                  uiOutput("run_history_table")
                              )
             )
    ),
    
    tabPanel("Monte Carlo",
             sidebarLayout(
               sidebarPanel(
                 numericInput("mc_runs","Monte Carlo runs",value=200,min=10,max=1000),
                 p("0.05 s recommended for speed/accuracy balance"),
                 sliderInput("chute_delay_std_dev",  "Ejection delay sd (s)",   value=1,  min=0.1,max=5),
                 sliderInput("launch_angle_std_dev", "Launch angle sd (deg)",   value=5,  min=0.1,max=10),
                 sliderInput("wind_speed_std_dev",   "Wind speed sd (%)",       value=5,  min=1,  max=99),
                 sliderInput("wind_dir_std_dev",     "Wind direction sd (deg)", value=30, min=1,  max=180),
                 sliderInput("dry_mass_std_dev",     "Dry mass sd (%)",         value=2,  min=0,  max=10),
                 sliderInput("cd_std_dev",           "Drag coefficient sd (%)", value=10, min=0,  max=30),
                 sliderInput("prop_mass_std_dev",    "Propellant mass sd (%)",  value=2,  min=0,  max=10),
                 sliderInput("montecarlo_precision", "Integration interval (s)",value=0.05,min=0.005,max=0.2)
               ),
               mainPanel(
                 leafletOutput("map",height=600),
                 p("Click map to set launch position"),
                 actionButton("run_mc","> Run Monte Carlo",class="btn-warning"),
                 conditionalPanel("output.has_mc",
                                  div(class="card",style="padding:14px;margin-top:12px;",
                                      verbatimTextOutput("landing_pct"))
                 )
               )
             )
    ),
    
    tabPanel("About",
             div(style="max-width:700px;margin:40px auto;",
                 div(class="card",style="padding:32px;",
                     h3(style="color:var(--c2);font-weight:800;letter-spacing:1px;margin-bottom:4px;","RRRocket 3D"),
                     p(style="color:var(--c4);font-size:0.8rem;letter-spacing:2px;text-transform:uppercase;margin-bottom:24px;",
                       "Model rocket flight simulator"),
                     tags$hr(style="border-color:#fc913a33;margin-bottom:24px;"),
                     h6(style="color:var(--c3);text-transform:uppercase;letter-spacing:1px;font-size:0.72rem;","What it does"),
                     p(style="color:var(--text);font-size:0.85rem;line-height:1.8;",
                       "RRRocket 3D simulates low-power model rocket flights in three dimensions,
         accounting for aerodynamic drag, motor thrust curves, wind weathercocking,
         parachute descent, and atmospheric density variation with altitude.
         A Monte Carlo engine propagates uncertainty in wind, ejection delay,
         and build tolerances to produce a probabilistic landing footprint."),
                     br(),
                     h6(style="color:var(--c3);text-transform:uppercase;letter-spacing:1px;font-size:0.72rem;","Physics"),
                     p(style="color:var(--text);font-size:0.85rem;line-height:1.8;",
                       "Aerodynamics use the extended Barrowman equations with a body-lift correction term.
         Drag is corrected for Mach number using a Prandtl-Glauert factor below Mach 0.8
         and a transonic ramp above it. Atmospheric density follows the
         International Standard Atmosphere model. Wind turbulence is modeled with
         an Ornstein-Uhlenbeck process driven by Weibull-distributed gust amplitudes."),
                     br(),
                     h6(style="color:var(--c3);text-transform:uppercase;letter-spacing:1px;font-size:0.72rem;","Monte Carlo"),
                     p(style="color:var(--text);font-size:0.85rem;line-height:1.8;",
                       "Each Monte Carlo run perturbs wind speed, wind direction, ejection delay,
         drag coefficient, dry mass, and propellant mass by user-specified standard
         deviations. Landing points are plotted on satellite imagery and scored
         against a user-drawn safe zone polygon. The 95th percentile drift radius
         is the recommended metric for NAR/Tripoli range safety submissions."),
                     br(),
                     h6(style="color:var(--c3);text-transform:uppercase;letter-spacing:1px;font-size:0.72rem;","Engine data"),
                     p(style="color:var(--text);font-size:0.85rem;line-height:1.8;",
                       "Thrust curves use the standard RASP .eng file format.
         Built-in curves are included for common Estes A-D motors. Upload any .eng file from",
                       tags$a(href="https://www.thrustcurve.org",target="_blank",style="color:var(--c2);","thrustcurve.org"),
                       "for higher-power motors."),
                     br(),
                     tags$hr(style="border-color:#fc913a33;margin:8px 0 20px;"),
                     div(style="display:flex;gap:12px;",
                         tags$a(href="https://github.com/tatecommission/rrrocket",target="_blank",class="btn btn-default","GitHub"),
                         tags$a(href="https://www.thrustcurve.org",target="_blank",class="btn btn-default","ThrustCurve.org")
                     )
                 )
             )
    ),
    
    nav_spacer(),
    nav_item(
      tags$a(
        href="https://github.com/tatecommission/rrrocket",
        target="_blank",
        style=paste0(
          "font-size:0.65rem;font-weight:700;letter-spacing:1px;",
          "text-transform:uppercase;color:#eae374;text-decoration:none;",
          "border:1px solid #fc913a55;border-radius:5px;padding:4px 10px;",
          "background:rgba(252,145,58,0.1);",
          "box-shadow:0 1px 0 rgba(255,255,255,0.08) inset,0 2px 5px rgba(0,0,0,0.4);",
          "display:inline-block;"
        ),
        "GitHub ↗"
      )
    )
  )
)

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

server <- function(input, output, session) {
  use_metric <- reactiveVal(TRUE)
  observeEvent(input$units, { use_metric(input$units == "metric") }, ignoreInit=TRUE)
  
  run_history   <- reactiveVal(list())
  results_store <- reactiveVal(NULL)
  mc_store      <- reactiveVal(NULL)
  landing_pct_val <- reactiveVal(NULL)
  
  output$has_results <- reactive({ !is.null(results_store()) })
  outputOptions(output, "has_results", suspendWhenHidden=FALSE)
  output$has_mc <- reactive({ !is.null(mc_store()) })
  outputOptions(output, "has_mc", suspendWhenHidden=FALSE)
  
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
  
  motor_data <- reactive({
    if (!is.null(input$motor_file)) {
      tryCatch(parse_thrust_input(input$motor_file, NULL), error=function(e) NULL)
    } else {
      ec <- input$engine_choice
      if (is.null(ec) || ec == "") return(NULL)
      tryCatch(parse_thrust_input(NULL, ec), error=function(e) NULL)
    }
  })
  
  aero_reactive <- reactive({
    if (!all(sapply(list(input$nose_type, input$fin_count), isTruthy))) return(NULL)
    if (!all(sapply(list(si$nose_length(), si$body_length(), si$diameter(),
                         si$fin_root(), si$fin_tip(), si$fin_span(),
                         si$fin_sweep(), si$cg_measured()),
                    function(x) isTruthy(x) && x > 0))) return(NULL)
    tryCatch(
      compute_aero(input$nose_type,
                   si$nose_length(), si$body_length(), si$diameter(), input$fin_count,
                   si$fin_root(), si$fin_tip(), si$fin_span(), si$fin_sweep(), si$cg_measured()),
      error=function(e) NULL)
  })
  
  aero_loaded <- reactive({
    aero <- aero_reactive(); if (is.null(aero)) return(NULL)
    td   <- motor_data()
    if (is.null(td)) return(c(aero, list(stability_margin_loaded=aero$stability_margin,
                                         cg_loaded=si$cg_measured())))
    cg_motor  <- si$nose_length() + si$body_length() - td$motor_length_m/2
    cg_loaded <- (si$dry_mass()*si$cg_measured() + td$prop_mass*cg_motor) /
      (si$dry_mass() + td$prop_mass)
    c(aero, list(stability_margin_loaded=(aero$CP-cg_loaded)/si$diameter(),
                 cg_loaded=cg_loaded))
  })
  
  output$stability_indicator <- renderUI({
    al <- aero_loaded(); if (is.null(al)) return(NULL)
    sm  <- al$stability_margin_loaded
    fmt <- function(m) sprintf("%.1f mm  /  %.2f in", m*1000, m*39.3701)
    cls <- if (sm<0.5) "stab-red" else if (sm<1.0) "stab-yellow" else if (sm<=3.0) "stab-green" else "stab-yellow"
    lbl <- if (sm<0.5) "UNSTABLE" else if (sm<1.0) "MARGINAL" else if (sm<=3.0) "STABLE" else "OVERSTABLE"
    hint <- if (sm<0.5) "Unstable. Move CG forward or increase fin size." else
      if (sm>3.0) "Overstable. Risk of weathercocking." else ""
    div(class=paste("stab-box", cls),
        tags$b(sprintf("%s — %.2f cal", lbl, sm)), tags$br(),
        sprintf("CP: %s", fmt(al$CP)), tags$br(),
        sprintf("CG: %s", fmt(al$cg_loaded)), tags$br(),
        if (!is.null(motor_data()))
          sprintf("Stability at burnout: %.2f cal", al$stability_margin)
        else tags$span(style="color:var(--dim);", "Load an engine to see loaded CG"),
        if (nchar(hint)>0) tagList(tags$br(), tags$span(hint)) else NULL)
  })
  
  output$thrust_curve_plot <- renderPlot({
    tc <- motor_data()$thrust_curve
    if (is.null(tc)||nrow(tc)==0) return(
      ggplot()+annotate("text",x=0.5,y=0.5,label="Select or upload an engine",
                        color="#6b7280",size=4)+theme_plot()+
        theme(axis.text=element_blank(),axis.title=element_blank(),panel.grid=element_blank()))
    tf <- approxfun(tc$time,tc$thrust,yleft=0,yright=0)
    ti <- integrate(tf,min(tc$time),max(tc$time))$value; mt <- max(tc$thrust)
    tc2 <- tc; if (!use_metric()) tc2$thrust <- tc2$thrust*N_to_lbf
    ylab  <- if (use_metric()) "thrust (N)" else "thrust (lbf)"
    t_ann <- if (use_metric()) sprintf("Total: %.2f Ns",ti) else sprintf("Total: %.2f lbf*s",ti*N_to_lbf)
    p_ann <- if (use_metric()) sprintf("Peak:  %.2f N",mt)  else sprintf("Peak:  %.2f lbf",mt*N_to_lbf)
    ggplot(tc2,aes(time,thrust))+geom_area(fill="#1a56db",alpha=0.08)+
      geom_line(color="#1a56db",linewidth=1)+geom_hline(yintercept=0,color="#e8eaf0")+
      annotate("text",x=max(tc2$time),y=max(tc2$thrust),label=t_ann,hjust=1,vjust=1.3,size=3.5,fontface="bold",color="#111928")+
      annotate("text",x=max(tc2$time),y=max(tc2$thrust)*0.87,label=p_ann,hjust=1,vjust=1.3,size=3.5,color="#6b7280")+
      scale_y_continuous(limits=c(0,max(tc2$thrust)*1.18))+
      labs(x="time (s)",y=ylab,title="Engine thrust vs. time")+theme_plot()
  })
  
  # ── SIMULATE ──────────────────────────────────────────────────────────────
  observeEvent(input$run, {
    aero   <- aero_reactive()
    parsed <- motor_data()
    if (is.null(aero))   { showNotification("Complete rocket geometry on Setup tab", type="warning"); return() }
    if (is.null(parsed)) { showNotification("Select an engine on Setup tab",         type="warning"); return() }
    sim <- tryCatch(
      flight_simulation_3d(
        parsed$thrust_curve, parsed$prop_mass, si$dry_mass(),
        si$diameter(), max(si$parachute_diam(),0.05), input$parachute_delay,
        aero$Cd, aero$CNa_total, aero$CP, input$precision,
        si$wind_speed(), input$wind_dir,
        si$cg_measured(), si$nose_length(), si$body_length(),
        parsed$motor_length_m, si$rail_length(),
        input$launch_bearing, input$launch_angle,
        landing_only=FALSE),
      error=function(e) { showNotification(paste("Sim error:",e$message),type="error"); NULL })
    if (is.null(sim)) return()
    res  <- list(sim=sim, aero=aero,
                 label=paste0("Run ",length(run_history())+1),
                 motor=if(isTruthy(input$engine_choice)&&input$engine_choice!="") input$engine_choice else "custom")
    hist <- run_history(); hist[[length(hist)+1]] <- res; run_history(hist)
    results_store(res)
  })
  
  # ── MONTE CARLO ───────────────────────────────────────────────────────────
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
            parsed$prop_mass * rnorm(1,1,0.01*input$prop_mass_std_dev),
            si$dry_mass()    * rnorm(1,1,0.01*input$dry_mass_std_dev),
            si$diameter(), max(si$parachute_diam(),0.05),
            max(0, rnorm(1,input$parachute_delay,input$chute_delay_std_dev)),
            aero$Cd * rnorm(1,1,0.01*input$cd_std_dev),
            aero$CNa_total, aero$CP, input$montecarlo_precision,
            max(0, rnorm(1,si$wind_speed(),0.01*si$wind_speed()*input$wind_speed_std_dev)),
            rnorm(1,input$wind_dir,input$wind_dir_std_dev),
            si$cg_measured(), si$nose_length(), si$body_length(),
            parsed$motor_length_m, si$rail_length(),
            input$launch_bearing,
            max(0, rnorm(1,input$launch_angle,input$launch_angle_std_dev)),
            landing_only=TRUE),
          error=function(e) NULL)
        landings[[i]] <- if (!is.null(sim)) data.frame(x=sim$x,y=sim$y) else data.frame(x=0,y=0)
      }
    })
    mc_store(do.call(rbind, landings))
  })
  
  # ── OUTPUTS ───────────────────────────────────────────────────────────────
  output$run_history_table <- renderUI({
    hist <- run_history(); if (length(hist)==0) return(p("No runs yet."))
    sc <- if(use_metric()) 1 else m_to_ft; u <- if(use_metric()) "m" else "ft"
    rows <- lapply(rev(seq_along(hist)), function(i) {
      r <- hist[[i]]; s <- r$sim
      tags$tr(tags$td(r$label), tags$td(r$motor),
              tags$td(sprintf("%.0f %s",max(s$altitude)*sc,u)),
              tags$td(sprintf("%.1f s",max(s$time))),
              tags$td(sprintf("%.2f cal",r$aero$stability_margin)))
    })
    tags$table(class="run-table",
               tags$thead(tags$tr(tags$th("Run"),tags$th("Motor"),
                                  tags$th("Apogee"),tags$th("Time"),tags$th("Stab."))),
               tags$tbody(rows))
  })
  
  output$summary <- renderPrint({
    res <- results_store(); req(!is.null(res))
    r <- res$sim; ae <- res$aero
    sc <- if(use_metric()) 1 else m_to_ft
    u  <- if(use_metric()) "m" else "ft"
    us <- if(use_metric()) "m/s" else "ft/s"
    cat(sprintf("apogee             %d %s\n",       round(max(r$altitude)*sc), u))
    cat(sprintf("max velocity       %.1f %s\n",     max(r$velocity)*sc, us))
    rail_v <- attr(r,"rail_exit_ms")
    if (!is.null(rail_v)&&!is.na(rail_v)) {
      flag <- if(rail_v*sc < 49) " *** LOW" else ""
      cat(sprintf("rail exit speed    %.1f %s%s\n", rail_v*sc, us, flag))
    }
    cat(sprintf("max Mach           %.3f\n",   max(r$mach)))
    cat(sprintf("time to apogee     %.2f s\n", r$time[which.max(r$altitude)]))
    cat(sprintf("total flight time  %.2f s\n", max(r$time)))
    cat(sprintf("stability (loaded) %.2f cal\n",round(min(r$stability_margin),2)))
    cat(sprintf("Cd                 %.4f\n",   ae$Cd))
  })
  
  output$altitude_plot <- renderPlot({
    res <- results_store(); req(!is.null(res)); r <- res$sim
    alt  <- if(use_metric()) r$altitude else r$altitude*m_to_ft
    ylab <- if(use_metric()) "altitude (m)" else "altitude (ft)"
    ggplot(data.frame(t=r$time,alt=alt),aes(t,alt))+
      geom_area(fill="#fc913a",alpha=0.15)+geom_line(color="#fc913a",linewidth=1)+
      geom_hline(yintercept=0,color="#e8eaf0")+
      labs(x="time (s)",y=ylab,title="Altitude")+theme_plot()
  })
  
  output$velocity_plot <- renderPlot({
    res <- results_store(); req(!is.null(res)); r <- res$sim
    vz   <- if(use_metric()) r$vz else r$vz*m_to_ft
    ylab <- if(use_metric()) "vertical velocity (m/s)" else "vertical velocity (ft/s)"
    ggplot(data.frame(t=r$time,vz=vz),aes(t,vz))+
      geom_line(color="#f9d62e",linewidth=1)+
      geom_hline(yintercept=0,color="#fc913a44",linetype="dashed")+
      labs(x="time (s)",y=ylab,title="Vertical velocity")+theme_plot()
  })
  
  output$track_3d <- renderPlotly({
    res <- results_store(); req(!is.null(res)); r <- res$sim
    sc <- if(use_metric()) 1 else m_to_ft
    xl <- if(use_metric()) "East (m)"     else "East (ft)"
    yl <- if(use_metric()) "North (m)"    else "North (ft)"
    zl <- if(use_metric()) "Altitude (m)" else "Altitude (ft)"
    plot_ly(r,x=~x*sc,y=~y*sc,z=~altitude*sc,type="scatter3d",mode="lines",
            line=list(color=~altitude*sc,colorscale=list(c(0,"#1a56db"),c(1,"#60a5fa")),width=3))|>
      add_trace(x=tail(r$x,1)*sc,y=tail(r$y,1)*sc,z=0,type="scatter3d",mode="markers",
                marker=list(color="#e02424",size=5),name="landing")|>
      layout(paper_bgcolor="#ffffff",font=list(color="#111928"),
             scene=list(bgcolor="#f8f9fb",
                        xaxis=list(title=xl,gridcolor="#d0d5de"),
                        yaxis=list(title=yl,gridcolor="#d0d5de"),
                        zaxis=list(title=zl,gridcolor="#d0d5de")))
  })
  
  output$landing_pct <- renderPrint({
    lc <- mc_store(); req(!is.null(lc))
    sc <- if(use_metric()) 1 else m_to_ft; u <- if(use_metric()) "m" else "ft"
    drift <- sqrt(lc$x^2+lc$y^2)*sc
    cat(sprintf("95th pct distance from pad: %.0f %s\n", quantile(drift,0.95), u))
    cat(sprintf("Max distance from pad:       %.0f %s\n", max(drift), u))
    if (!is.null(landing_pct_val()))
      cat(sprintf("%.1f%% of rockets land in the polygon\n", landing_pct_val()))
  })
  
  # ── MAP ───────────────────────────────────────────────────────────────────
  launch_point <- reactiveVal(list(lat=38.89, lng=-77.03))
  observeEvent(input$map_click, {
    launch_point(list(lat=input$map_click$lat, lng=input$map_click$lng))
  })
  
  drawn_polygon <- reactiveVal(NULL)
  observeEvent(input$map_draw_new_feature, { drawn_polygon(input$map_draw_new_feature) })
  
  rocket_icon <- makeIcon(
    iconUrl="https://raw.githubusercontent.com/pointhi/leaflet-color-markers/master/img/marker-icon-2x-red.png",
    iconWidth=25,iconHeight=41,iconAnchorX=12,iconAnchorY=41)
  
  output$map <- renderLeaflet({
    leaflet()|>addProviderTiles("Esri.WorldImagery")|>
      setView(lat=36.9052,lng=-81.0768,zoom=5)|>
      addDrawToolbar(polylineOptions=FALSE,circleOptions=drawCircleOptions(),
                     markerOptions=FALSE,circleMarkerOptions=FALSE,
                     rectangleOptions=drawRectangleOptions(),
                     polygonOptions=drawPolygonOptions(),editOptions=FALSE)
  })
  
  observe({
    lp <- launch_point()
    leafletProxy("map")|>clearGroup("launch")|>
      addMarkers(lng=lp$lng,lat=lp$lat,icon=rocket_icon,label="Launch pad",group="launch")
  })
  
  observeEvent(list(mc_store(), drawn_polygon()), {
    lc <- mc_store(); poly <- drawn_polygon()
    if (is.null(lc)||is.null(poly)) return()
    lp <- launch_point(); lat0 <- lp$lat; lng0 <- lp$lng
    lc$lat <- lat0+(lc$y/111320)
    lc$lng <- lng0+(lc$x/(111320*cos(lat0*pi/180)))
    poly_coords <- poly$geometry$coordinates[[1]]
    poly_mat <- do.call(rbind,lapply(poly_coords,function(p) c(p[[1]],p[[2]])))
    if (!identical(poly_mat[1,],poly_mat[nrow(poly_mat),])) poly_mat <- rbind(poly_mat,poly_mat[1,])
    pip <- function(px,py,pm) sapply(seq_along(px),function(k){
      x<-px[k];y<-py[k];n<-nrow(pm);j<-n;inside<-FALSE
      for(i in 1:n){xi<-pm[i,1];yi<-pm[i,2];xj<-pm[j,1];yj<-pm[j,2]
      if(((yi>y)!=(yj>y))&&(x<(xj-xi)*(y-yi)/(yj-yi)+xi)) inside<-!inside;j<-i}
      inside})
    inside <- pip(lc$lng,lc$lat,poly_mat)
    pct    <- mean(inside)*100
    landing_pct_val(pct)
    lc$color <- ifelse(inside,"#0e9f6e","#e02424")
    icon <- makeIcon(iconUrl="https://raw.githubusercontent.com/pointhi/leaflet-color-markers/master/img/marker-icon-2x-red.png",
                     iconWidth=25,iconHeight=41,iconAnchorX=12,iconAnchorY=41)
    leafletProxy("map")|>clearMarkers()|>clearGroup("launch")|>
      addCircleMarkers(lng=lc$lng,lat=lc$lat,radius=4,color=lc$color,fillOpacity=0.7,stroke=FALSE)|>
      addMarkers(lng=lng0,lat=lat0,icon=icon,label="Launch pad",group="launch")|>
      addPopups(lng=lng0,lat=lat0,popup=sprintf("%.1f%% land inside safe zone",pct))
  })
}

shinyApp(ui=ui, server=server)