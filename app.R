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
    pm    <- as.numeric(hdr[5])
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
  --bg:     #1f1f1f;
  --panel:  #2a2a2a;
  --border: #3a3a3a;
  --accent: #5ea0ff;
  --accent2:#34d399;
  --warn:   #f87171;
  --yellow: #fbbf24;
  --text:   #ffffff;
  --dim:    #d9d9d9;
  --sans:   'Lexend', sans-serif;
}

*,*::before,*::after {
  font-family: 'Lexend', sans-serif !important;
}

html,body{
  background:var(--bg)!important;
  color:var(--text)!important;
  font-family:var(--sans)!important;
  font-weight:400!important;
}

.navbar{
  background:var(--panel)!important;
  border-bottom:1px solid var(--border)!important;
  padding:0 20px!important;
  box-shadow:0 1px 3px rgba(0,0,0,.35)!important;
}

.navbar-brand{
  font-size:1rem!important;
  font-weight:700!important;
  color:var(--text)!important;
  letter-spacing:.5px!important;
  font-family:var(--sans)!important;
}

.nav-link{
  font-size:0.78rem!important;
  color:var(--dim)!important;
  padding:14px 16px!important;
  border-bottom:2px solid transparent!important;
  text-transform:uppercase!important;
  letter-spacing:.8px!important;
  font-weight:600!important;
  font-family:var(--sans)!important;
}

.nav-link:hover{
  color:var(--text)!important;
}

.nav-link.active{
  color:var(--accent)!important;
  border-bottom-color:var(--accent)!important;
  background:transparent!important;
}

.tab-content{
  padding:20px 24px;
}

label,.form-label{
  font-size:0.72rem!important;
  color:var(--text)!important;
  text-transform:uppercase!important;
  letter-spacing:.05em!important;
  font-weight:600!important;
  font-family:var(--sans)!important;
}

.form-control,.form-select,input[type=number],select{
  background:var(--panel)!important;
  border:1px solid var(--border)!important;
  color:var(--text)!important;
  font-family:var(--sans)!important;
  font-size:0.85rem!important;
  font-weight:500!important;
  border-radius:4px!important;
  padding:6px 10px!important;
}

.form-control:focus,input:focus,select:focus{
  border-color:var(--accent)!important;
  box-shadow:0 0 0 3px rgba(94,160,255,.15)!important;
  outline:none!important;
}

select option{
  background:var(--panel);
  color:var(--text);
}

.irs--shiny .irs-bar{
  background:var(--accent)!important;
  border-color:var(--accent)!important;
}

.irs--shiny .irs-handle{
  border-color:var(--accent)!important;
  background:#fff!important;
}

.irs--shiny .irs-single{
  background:var(--accent)!important;
  color:#fff!important;
  font-family:var(--sans)!important;
  font-size:0.68rem!important;
  font-weight:600!important;
}

.irs--shiny .irs-line{
  background:var(--border)!important;
}

.irs--shiny .irs-grid-text,
.irs--shiny .irs-min,
.irs--shiny .irs-max{
  color:#ffffff!important;
  font-family:var(--sans)!important;
  font-weight:500!important;
  background:transparent!important;
}

.irs--shiny .irs-grid-pol.small,
.irs--shiny .irs-grid-pol{
  background:#ffffff!important;
}

.btn{
  font-size:0.75rem!important;
  text-transform:uppercase!important;
  letter-spacing:.8px!important;
  border-radius:4px!important;
  padding:7px 18px!important;
  font-weight:700!important;
  font-family:var(--sans)!important;
  transition:all .15s!important;
}

.btn-primary{
  background:var(--accent)!important;
  border-color:var(--accent)!important;
  color:#fff!important;
}

.btn-primary:hover{
  background:#3f83f8!important;
  border-color:#3f83f8!important;
}

.btn-warning{
  background:transparent!important;
  border:1px solid var(--accent)!important;
  color:var(--accent)!important;
}

.btn-warning:hover{
  background:rgba(94,160,255,.08)!important;
}

.btn-default,.btn-secondary{
  background:transparent!important;
  border:1px solid var(--border)!important;
  color:var(--dim)!important;
}

p{
  font-size:0.72rem;
  color:var(--dim);
  margin:2px 0 10px;
  font-family:var(--sans)!important;
  font-weight:400!important;
}

.shiny-plot-output{
  background:var(--panel)!important;
  border:1px solid var(--border);
  border-radius:4px;
}

pre,.shiny-verbatim-output{
  background:#181818!important;
  color:#f5f5f5!important;
  border:1px solid var(--border)!important;
  border-radius:4px!important;
  font-family:var(--sans)!important;
  font-size:0.8rem!important;
  padding:14px!important;
  line-height:1.8!important;
}

.leaflet-container{
  border-radius:4px;
  border:1px solid var(--border);
}

.well{
  background:var(--panel)!important;
  border:1px solid var(--border)!important;
  border-radius:4px!important;
  box-shadow:none!important;
}

.nav-pills .nav-link{
  font-size:0.7rem!important;
  color:var(--dim)!important;
  border-radius:4px!important;
  padding:5px 12px!important;
  text-transform:uppercase!important;
  letter-spacing:.05em!important;
  background:transparent!important;
  border:1px solid var(--border)!important;
  margin:2px!important;
  font-weight:600!important;
  font-family:var(--sans)!important;
}

.nav-pills .nav-link.active{
  background:var(--accent)!important;
  border-color:var(--accent)!important;
  color:#fff!important;
}

.card{
  background:var(--panel)!important;
  border:1px solid var(--border)!important;
  border-radius:4px!important;
  box-shadow:0 1px 3px rgba(0,0,0,.25)!important;
}

.card-body{
  padding:16px!important;
}

.home-wrap{
  display:flex;
  flex-direction:column;
  justify-content:center;
  align-items:center;
  height:70vh;
  text-align:center;
}

.title-wrap{
  display:flex;
  align-items:center;
  justify-content:center;
}

.home-title{
  font-size:3.4rem;
  font-weight:700;
  color:var(--text);
  letter-spacing:1px;
  margin:0;
  font-family:var(--sans)!important;
  text-rendering:geometricPrecision;
}

.cursor{
  display:inline-block;
  width:3px;
  height:3.4rem;
  background:white;
  margin-left:4px;
  vertical-align:middle;
  animation:blink 0.7s step-end infinite;
}

@keyframes blink { 50% { opacity:0; } }

.home-sub{
  font-size:0.85rem;
  color:var(--dim);
  letter-spacing:2px;
  margin-top:8px;
  text-transform:uppercase;
  font-family:var(--sans)!important;
  font-weight:500!important;
}

.stab-box{
  border:1px solid var(--border);
  border-radius:4px;
  padding:10px 14px;
  margin-top:10px;
  font-family:var(--sans)!important;
  font-size:0.78rem;
  line-height:1.9;
  background:var(--panel);
  font-weight:500!important;
}

.stab-green{
  color:#86efac;
  border-color:#22c55e!important;
  background:#132218!important;
}

.stab-yellow{
  color:#fde68a;
  border-color:#f59e0b!important;
  background:#2b2111!important;
}

.stab-red{
  color:#fca5a5;
  border-color:#ef4444!important;
  background:#2a1515!important;
}

.run-table{
  width:100%;
  border-collapse:collapse;
  font-size:0.75rem;
  margin-top:8px;
  font-family:var(--sans)!important;
}

.run-table th{
  color:var(--dim);
  text-transform:uppercase;
  letter-spacing:.05em;
  padding:5px 8px;
  border-bottom:2px solid var(--border);
  text-align:left;
  font-weight:700;
}

.run-table td{
  color:var(--text);
  padding:5px 8px;
  border-bottom:1px solid var(--border);
  font-weight:500;
}

.run-table tr:hover td{
  background:#303030;
}
"

unit_toggle_ui <- function(id) {
  div(style="margin-bottom:14px;",
      radioButtons(id,label=NULL,choices=c("Metric"="metric","Imperial"="imperial"),
                   selected="metric",inline=TRUE))
}

labeled_num <- function(input_id, label_id, default, min=0) {
  div(style="margin-bottom:8px;",
      uiOutput(label_id,inline=FALSE),
      numericInput(input_id,label=NULL,value=default,min=min,width="100%"))
}

ui <- tagList(
  tags$head(tags$style(HTML(css))),
  navbarPage(
    title="RRRocket 3D",
    theme=bs_theme(version=5,bg="#f5f6f8",fg="#111928",primary="#1a56db"),
    
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
                 p(class="home-sub","Model rocket flight simulator"),
                 br(), unit_toggle_ui("units")
             )
    ),
    
    tabPanel("Setup",
             fluidRow(
               column(5,
                      unit_toggle_ui("units2"),
                      navset_card_pill(
                        nav_panel("Rocket",
                                  labeled_num("dry_mass_val","lbl_dry_mass",90),
                                  labeled_num("diameter","lbl_diameter",24),
                                  labeled_num("body_length","lbl_body_length",300),
                                  labeled_num("cg_measured","lbl_cg_measured",220)),
                        nav_panel("Nosecone",
                                  selectInput("nose_type","Nosecone type",choices=c("ogive","conical","parabolic")),
                                  labeled_num("nose_length","lbl_nose_length",70)),
                        nav_panel("Fins",
                                  numericInput("fin_count","Number of fins",value=3),
                                  labeled_num("fin_root","lbl_fin_root",50),
                                  labeled_num("fin_tip","lbl_fin_tip",25),
                                  labeled_num("fin_span","lbl_fin_span",30),
                                  labeled_num("fin_sweep","lbl_fin_sweep",20)),
                        nav_panel("Chute",
                                  labeled_num("parachute_diameter","lbl_chute_diameter",305)),
                        nav_panel("Launch",
                                  labeled_num("rail_length","lbl_rail_length",0.9),
                                  labeled_num("wind_speed_val","lbl_wind_speed",3),
                                  sliderInput("wind_dir","Wind from (deg CW from N)",min=0,max=360,value=270),
                                  numericInput("launch_angle","Launch angle from vertical (deg)",value=0,min=0,max=30),
                                  sliderInput("launch_bearing","Launch bearing (deg CW from N)",value=0,min=0,max=360)),
                        nav_panel("Engine",
                                  fileInput("motor_file",NULL,accept=".eng",buttonLabel="Upload .eng"),
                                  numericInput("parachute_delay","Ejection delay (s)",value=4),
                                  selectInput("engine_choice","Or choose an engine:",
                                              choices=c("Select engine..."="","A8","A10","B4","B6","C6","C11","D12","E12","E16",
                                              "G40","L2350"),
                                              selected="B6",size=5,selectize=FALSE))),
                      uiOutput("stability_indicator")),
               column(7,plotOutput("thrust_curve_plot",height="420px")))),
    
    tabPanel("Simulation",
             sliderInput("precision","Integration interval (s)",value=0.01,min=0.001,max=0.1),
             p("0.01 s recommended for accuracy"),
             actionButton("run","Simulate",class="btn-primary"),
             br(),br(),
             h5(style="color:var(--dim);font-family:var(--sans);font-size:0.8rem;text-transform:uppercase;letter-spacing:1px;","Run history"),
             uiOutput("run_history_table")),
    
    tabPanel("Results",
             plotOutput("altitude_plot"),
             plotOutput("velocity_plot"),
             verbatimTextOutput("summary")),
    
    tabPanel("3D Track",
             plotlyOutput("track_3d",width="100%",height="780px"),
             verbatimTextOutput("landing_summary")),
    
    tabPanel("Monte Carlo",
             sidebarLayout(
               sidebarPanel(
                 numericInput("mc_runs","Monte Carlo runs",value=200,min=10,max=1000),
                 p("0.05 s recommended for speed/accuracy balance"),
                 sliderInput("chute_delay_std_dev",  "Ejection delay sd (s)",      value=1,  min=0.1,max=5),
                 sliderInput("launch_angle_std_dev", "Launch angle sd (deg)",      value=5,  min=0.1,max=10),
                 sliderInput("wind_speed_std_dev",   "Wind speed sd (%)",          value=5,  min=1,  max=99),
                 sliderInput("wind_dir_std_dev",     "Wind direction sd (deg)",    value=30, min=1,  max=180),
                 sliderInput("dry_mass_std_dev",     "Dry mass sd (%)",            value=2,  min=0,  max=10),
                 sliderInput("cd_std_dev",           "Drag coefficient sd (%)",    value=10, min=0,  max=30),
                 sliderInput("prop_mass_std_dev",    "Propellant mass sd (%)",     value=2,  min=0,  max=10)),
               mainPanel(
                 leafletOutput("map",height=600),
                 p("Click map to set launch position"),
                 actionButton("run_mc","Run Monte Carlo",class="btn-warning"),
                 verbatimTextOutput("landing_pct"))))
  )
)

theme_plot <- function() {
  theme_minimal(base_size=11)+theme(
    plot.background  = element_rect(fill="#ffffff",color=NA),
    panel.background = element_rect(fill="#ffffff",color=NA),
    panel.grid.major = element_line(color="#e8eaf0",linewidth=0.4),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color="#d0d5de",fill=NA,linewidth=0.5),
    axis.text        = element_text(color="#6b7280",size=8),
    axis.title       = element_text(color="#374151",size=9),
    plot.title       = element_text(color="#111928",size=11,face="bold"),
    plot.margin      = margin(8,12,8,8))
}

server <- function(input, output, session) {
  
  use_metric <- reactiveVal(TRUE)
  observeEvent(input$units,  { use_metric(input$units=="metric");  updateRadioButtons(session,"units2",selected=input$units)  },ignoreInit=TRUE)
  observeEvent(input$units2, { use_metric(input$units2=="metric"); updateRadioButtons(session,"units", selected=input$units2) },ignoreInit=TRUE)
  
  mk_lbl <- function(id, text, um, ui) {
    output[[id]] <- renderUI(tags$label(
      style="font-family:var(--sans);font-size:0.72rem;color:var(--dim);text-transform:uppercase;letter-spacing:.05em;",
      paste0(text," (",if(use_metric()) um else ui,")")))
  }
  mk_lbl("lbl_dry_mass",       "Dry mass",               "g",   "oz")
  mk_lbl("lbl_diameter",       "Body tube diameter",     "mm",  "in")
  mk_lbl("lbl_body_length",    "Body tube length",       "mm",  "in")
  mk_lbl("lbl_cg_measured",    "CG from nose tip (dry)", "mm",  "in")
  mk_lbl("lbl_nose_length",    "Nosecone length",        "mm",  "in")
  mk_lbl("lbl_fin_root",       "Root chord",             "mm",  "in")
  mk_lbl("lbl_fin_tip",        "Tip chord",              "mm",  "in")
  mk_lbl("lbl_fin_span",       "Semi-span",              "mm",  "in")
  mk_lbl("lbl_fin_sweep",      "Sweep length",           "mm",  "in")
  mk_lbl("lbl_chute_diameter", "Chute diameter",         "mm",  "in")
  mk_lbl("lbl_rail_length",    "Rail length",            "m",   "ft")
  mk_lbl("lbl_wind_speed",     "Wind speed",             "m/s", "mph")
  
  to_m   <- function(v) { if(!isTruthy(v)) return(0); if(use_metric()) v/1000 else v*0.0254 }
  to_m_r <- function(v) { if(!isTruthy(v)) return(0); if(use_metric()) v      else v*0.3048 }
  get_dry_mass_kg  <- function() { v<-input$dry_mass_val;       if(!isTruthy(v)) return(0.09); if(use_metric()) v/1000 else v*0.0283495 }
  get_chute_diam_m <- function() { v<-input$parachute_diameter; if(!isTruthy(v)) return(0.30); if(use_metric()) v/1000 else v*0.0254 }
  get_wind_ms      <- function() { v<-input$wind_speed_val;     if(!isTruthy(v)) return(3);    if(use_metric()) v      else v*0.44704 }
  
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
    needed <- list(input$nose_type,input$nose_length,input$diameter,input$body_length,
                   input$cg_measured,input$fin_count,input$fin_root,input$fin_tip,
                   input$fin_span,input$fin_sweep)
    if(!all(sapply(needed,isTruthy))) return(NULL)
    tryCatch(compute_aero(input$nose_type,
                          to_m(input$nose_length),to_m(input$body_length),
                          to_m(input$diameter),input$fin_count,
                          to_m(input$fin_root),to_m(input$fin_tip),
                          to_m(input$fin_span),to_m(input$fin_sweep),
                          to_m(input$cg_measured)),
             error=function(e) NULL)
  })
  
  aero_loaded <- reactive({
    aero <- aero_reactive(); if(is.null(aero)) return(NULL)
    td   <- motor_data()
    if(is.null(td)) return(c(aero,list(stability_margin_loaded=aero$stability_margin,
                                       cg_loaded=to_m(input$cg_measured))))
    prop_mass <- td$prop_mass; dry_mass <- get_dry_mass_kg()
    cg_dry    <- to_m(input$cg_measured)
    cg_motor  <- to_m(input$nose_length)+to_m(input$body_length)-td$motor_length_m/2
    cg_loaded <- (dry_mass*cg_dry+prop_mass*cg_motor)/(dry_mass+prop_mass)
    sm_loaded <- (aero$CP-cg_loaded)/to_m(input$diameter)
    c(aero,list(stability_margin_loaded=sm_loaded,cg_loaded=cg_loaded))
  })
  
  output$stability_indicator <- renderUI({
    al <- aero_loaded()
    if (is.null(al)) return(NULL)
    sm   <- al$stability_margin_loaded
    fmt  <- function(m) if (use_metric()) sprintf("%.0f mm", m*1000) else sprintf("%.2f in", m*39.37)
    cls  <- if (sm < 0.5) "stab-red" else if (sm < 1.0) "stab-yellow" else if (sm <= 3.0) "stab-green" else "stab-yellow"
    lbl  <- if (sm < 0.5) "UNSTABLE" else if (sm < 1.0) "MARGINAL" else if (sm <= 3.0) "STABLE" else "OVERSTABLE"
    hint <- if (sm < 0.5) "Unstable. Move CG forward or increase fin size." else
      if (sm > 3.0) "Overstable. Risk of weathercocking." else ""
    div(class = paste("stab-box", cls),
        tags$b(sprintf("%s — %.2f cal", lbl, sm)),
        tags$br(),
        sprintf("CP: %s   CG: %s", fmt(al$CP), fmt(al$cg_loaded)),
        tags$br(),
        if (!is.null(motor_data()))
          sprintf("Stability at burnout: %.2f cal", al$stability_margin)
        else
          tags$span(style="color:var(--dim);", "Load an engine to see loaded CG"),
        if (nchar(hint) > 0) tagList(tags$br(), tags$span(hint)) else NULL
    )
  })
  
  thrust_data <- reactive({ md<-motor_data(); if(is.null(md)) NULL else md$thrust_curve })
  
  output$thrust_curve_plot <- renderPlot({
    tc <- thrust_data()
    if(is.null(tc)||nrow(tc)==0) return(
      ggplot()+annotate("text",x=0.5,y=0.5,label="Select or upload an engine",
                        color="#6b7280",size=4)+theme_plot()+
        theme(axis.text=element_blank(),axis.title=element_blank(),panel.grid=element_blank()))
    tf <- approxfun(tc$time,tc$thrust,yleft=0,yright=0)
    ti <- integrate(tf,min(tc$time),max(tc$time))$value; mt <- max(tc$thrust)
    tc2 <- tc; if(!use_metric()) tc2$thrust <- tc2$thrust*N_to_lbf
    ylab  <- if(use_metric())"thrust (N)"else"thrust (lbf)"
    t_ann <- if(use_metric()) sprintf("Total: %.2f Ns",ti) else sprintf("Total: %.2f lbf*s",ti*N_to_lbf)
    p_ann <- if(use_metric()) sprintf("Peak:  %.2f N",mt)  else sprintf("Peak:  %.2f lbf",mt*N_to_lbf)
    ggplot(tc2,aes(time,thrust))+geom_area(fill="#1a56db",alpha=0.08)+
      geom_line(color="#1a56db",linewidth=1)+geom_hline(yintercept=0,color="#e8eaf0")+
      annotate("text",x=max(tc2$time),y=max(tc2$thrust),label=t_ann,hjust=1,vjust=1.3,
               size=3.5,fontface="bold",color="#111928")+
      annotate("text",x=max(tc2$time),y=max(tc2$thrust)*0.87,label=p_ann,hjust=1,vjust=1.3,
               size=3.5,color="#6b7280")+
      scale_y_continuous(limits=c(0,max(tc2$thrust)*1.18))+
      labs(x="time (s)",y=ylab,title="Engine thrust vs. time")+theme_plot()
  })
  
  run_history <- reactiveVal(list())
  
  results <- eventReactive(input$run, {
    aero   <- aero_reactive()
    parsed <- motor_data()
    validate(need(!is.null(aero),  "Complete rocket geometry on the Setup tab."))
    validate(need(!is.null(parsed),"Select an engine on the Setup tab."))
    tryCatch({
      sim <- flight_simulation_3d(
        parsed$thrust_curve,parsed$prop_mass,get_dry_mass_kg(),
        to_m(input$diameter),get_chute_diam_m(),input$parachute_delay,
        aero$Cd,aero$CNa_total,aero$CP,input$precision,
        get_wind_ms(),input$wind_dir,
        to_m(input$cg_measured),to_m(input$nose_length),to_m(input$body_length),
        parsed$motor_length_m,to_m_r(input$rail_length),
        input$launch_bearing,input$launch_angle,
        landing_only=FALSE)
      res <- list(sim=sim,aero=aero,
                  label=paste0("Run ",length(run_history())+1),
                  motor=if(isTruthy(input$engine_choice)&&input$engine_choice!="") input$engine_choice else "custom")
      hist <- run_history(); hist[[length(hist)+1]] <- res; run_history(hist)
      res
    },error=function(e){showNotification(paste("Error:",e$message),type="error",duration=8);NULL})
  })
  
  output$run_history_table <- renderUI({
    hist <- run_history()
    if(length(hist)==0) return(p("No runs yet."))
    rows <- lapply(rev(seq_along(hist)),function(i) {
      r <- hist[[i]]; s <- r$sim
      sc <- if(use_metric()) 1 else m_to_ft; u <- if(use_metric())"m"else"ft"
      tags$tr(tags$td(r$label),tags$td(r$motor),
              tags$td(sprintf("%.0f %s",max(s$altitude)*sc,u)),
              tags$td(sprintf("%.1f s",max(s$time))),
              tags$td(sprintf("%.2f cal",r$aero$stability_margin)))
    })
    tags$table(class="run-table",
               tags$thead(tags$tr(tags$th("Run"),tags$th("Motor"),tags$th("Apogee"),
                                  tags$th("Flight time"),tags$th("Stability"))),
               tags$tbody(rows))
  })
  
  output$altitude_plot <- renderPlot({
    req(results()); r <- results()$sim
    alt  <- if(use_metric()) r$altitude else r$altitude*m_to_ft
    ylab <- if(use_metric())"altitude (m)"else"altitude (ft)"
    ggplot(data.frame(t=r$time,alt=alt),aes(t,alt))+
      geom_area(fill="#1a56db",alpha=0.07)+geom_line(color="#1a56db",linewidth=1)+
      geom_hline(yintercept=0,color="#e8eaf0")+
      labs(x="time (s)",y=ylab,title="Altitude")+theme_plot()
  })
  
  output$velocity_plot <- renderPlot({
    req(results()); r <- results()$sim
    vz   <- if(use_metric()) r$vz else r$vz*m_to_ft
    ylab <- if(use_metric())"vertical velocity (m/s)"else"vertical velocity (ft/s)"
    ggplot(data.frame(t=r$time,vz=vz),aes(t,vz))+
      geom_line(color="#0e9f6e",linewidth=1)+
      geom_hline(yintercept=0,color="#d0d5de",linetype="dashed")+
      labs(x="time (s)",y=ylab,title="Vertical velocity")+theme_plot()
  })
  
  output$track_3d <- renderPlotly({
    req(results()); r <- results()$sim
    sc <- if(use_metric()) 1 else m_to_ft
    xl <- if(use_metric())"East (m)"else"East (ft)"
    yl <- if(use_metric())"North (m)"else"North (ft)"
    zl <- if(use_metric())"Altitude (m)"else"Altitude (ft)"
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
  
  output$summary <- renderPrint({
    req(results()); r <- results()$sim; ae <- results()$aero
    sc <- if(use_metric()) 1 else m_to_ft
    u  <- if(use_metric())"m"else"ft"; us <- if(use_metric())"m/s"else"ft/s"
    cat("apogee            ",round(max(r$altitude)*sc),u,"\n")
    cat("max velocity      ",round(max(r$velocity)*sc,1),us,"\n")
    rail_v <- attr(r,"rail_exit_ms")
    if(!is.null(rail_v)&&!is.na(rail_v)) {
      flag <- if(rail_v<15)" *** LOW - risk of instability"else""
      cat(sprintf("rail exit speed   %.1f %s%s\n",rail_v*sc,us,flag))
    }
    cat("max Mach number ",round(max(r$mach),3),"\n")
    cat("time at apogee ",round(r$time[which.max(r$altitude)],2),"s\n")
    cat("total flight time ",round(max(r$time),2),"s\n")
    cat("static stability when full ",round(ae$stability_margin,2),"calibers\n")
    cat("stability margin when loaded ",round(min(r$stability_margin),2),"cal\n")
    cat("Cd (subsonic) ",round(ae$Cd,4),"\n")
    cat("  nose",round(ae$Cd_nose,4)," fins",round(ae$Cd_fins,4),
        " body",round(ae$Cd_body,4)," below engine",round(ae$Cd_base,4),"\n")
  })
  
  output$landing_summary <- renderPrint({
    req(results()); r <- results()$sim
    sc <- if(use_metric()) 1 else m_to_ft; u <- if(use_metric())"m"else"ft"
    lx <- tail(r$x,1)*sc; ly <- tail(r$y,1)*sc
    cat(sprintf("East drift   %+.1f %s\n",lx,u))
    cat(sprintf("North drift  %+.1f %s\n",ly,u))
    cat(sprintf("Total drift  %.1f %s\n",sqrt(lx^2+ly^2),u))
  })
  
  landing_pct_val   <- reactiveVal(NULL)
  mc_landings_store <- reactiveVal(NULL)
  
  monte_carlo <- eventReactive(input$run_mc, {
    aero   <- aero_reactive()
    parsed <- motor_data()
    validate(need(!is.null(aero),   "Complete rocket geometry first."))
    validate(need(!is.null(parsed), "Select an engine first."))
    n        <- input$mc_runs
    landings <- vector("list", n)
    
    result <- tryCatch({
      withProgress(message="Monte Carlo", value=0, {
        for (i in seq_len(n)) {
          incProgress(1/n, detail=sprintf("Run %d / %d", i, n))
          dm_kg  <- get_dry_mass_kg() * rnorm(1, 1, 0.01*input$dry_mass_std_dev)
          pm_kg  <- parsed$prop_mass  * rnorm(1, 1, 0.01*input$prop_mass_std_dev)
          cd_var <- aero$Cd           * rnorm(1, 1, 0.01*input$cd_std_dev)
          sim <- tryCatch(
            flight_simulation_3d(
              parsed$thrust_curve, pm_kg, dm_kg,
              to_m(input$diameter), get_chute_diam_m(),
              max(0, rnorm(1, input$parachute_delay, input$chute_delay_std_dev)),
              cd_var, aero$CNa_total, aero$CP, input$precision,
              max(0, rnorm(1, get_wind_ms(), 0.01*get_wind_ms()*input$wind_speed_std_dev)),
              rnorm(1, input$wind_dir, input$wind_dir_std_dev),
              to_m(input$cg_measured), to_m(input$nose_length), to_m(input$body_length),
              parsed$motor_length_m, to_m_r(input$rail_length),
              input$launch_bearing,
              max(0, rnorm(1, input$launch_angle, input$launch_angle_std_dev)),
              landing_only=TRUE),  
            error=function(e) NULL)
          landings[[i]] <- if (!is.null(sim))
            data.frame(x=sim$x, y=sim$y) 
          else
            data.frame(x=0, y=0)
        }
      })
      landings
    },
    error = function(e) {
     
      completed <- Filter(Negate(is.null), landings)
      if (length(completed) > 0) {
        showNotification(sprintf("Stopped early — showing %d runs.", length(completed)), 
                         type="warning", duration=4)
        completed
      } else {
        NULL
      }
    })
    
    if (is.null(result)) return(NULL)
    lc <- do.call(rbind, result)
    mc_landings_store(lc)
    lc
  })
  
  eval_polygon <- function(lc, poly) {
    req(launch_point())
    lp <- launch_point(); lat0 <- lp$lat; lng0 <- lp$lng
    lc$lat <- lat0+(lc$y/111320)
    lc$lng <- lng0+(lc$x/(111320*cos(lat0*pi/180)))
    poly_coords <- poly$geometry$coordinates[[1]]
    poly_mat    <- do.call(rbind,lapply(poly_coords,function(p) c(p[[1]],p[[2]])))
    if(!identical(poly_mat[1,],poly_mat[nrow(poly_mat),])) poly_mat <- rbind(poly_mat,poly_mat[1,])
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
  }
  
  observeEvent(list(monte_carlo(),drawn_polygon()), {
    lc <- mc_landings_store(); poly <- drawn_polygon()
    if(is.null(lc)||is.null(poly)) return()
    eval_polygon(lc,poly)
  })
  
  output$landing_pct <- renderPrint({
    req(monte_carlo())
    lc    <- monte_carlo()
    sc    <- if(use_metric()) 1 else m_to_ft; u <- if(use_metric())"m"else"ft"
    drift <- sqrt(lc$x^2+lc$y^2)*sc
    cat(sprintf("95th pct distance from pad: %.0f %s\n",quantile(drift,0.95),u))
    cat(sprintf("Max distance from pad:       %.0f %s\n",max(drift),u))
    if(!is.null(landing_pct_val())) cat(sprintf("%.1f%% of rockets land in the polygon\n",landing_pct_val()))
  })
  
  launch_point <- reactiveVal(list(lat=38.89,lng=-77.03))
  observeEvent(input$map_click,{launch_point(list(lat=input$map_click$lat,lng=input$map_click$lng))})
  
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
  
  drawn_polygon <- reactiveVal(NULL)
  observeEvent(input$map_draw_new_feature,{drawn_polygon(input$map_draw_new_feature)})
}

shinyApp(ui=ui,server=server)