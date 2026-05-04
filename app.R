library(shiny)
library(bslib)
library(plotly)
library(ggplot2)
library(sf)
library(leaflet)
library(leaflet.extras)
library(shinylive)
library(httpuv)

g         <- 9.81
rho       <- 1.225
m_to_ft <- 3.28084


# DEFINE COMPUTE AERO

compute_aero <- function(nose_type, nose_length, body_length,
                         bt_diameter, fin_count, fin_root,
                         fin_tip, fin_span, fin_sweep, cg_measured) {
  
  bt_radius <- bt_diameter / 2
  Aref      <- pi * bt_radius^2
  
  # nose cone
  Xcp_nose <- switch(nose_type,
                     conical   = (2/3) * nose_length,
                     ogive     = 0.466 * nose_length,
                     parabolic = 0.5   * nose_length
  )
  
  half_angle <- atan(bt_radius / nose_length)
  
  Cd_nose <- switch(nose_type,
                    conical   = 0.8 * sin(half_angle)^2,
                    ogive     = 0.5 * sin(half_angle)^2,
                    parabolic = 0.3 * (bt_diameter / nose_length)^2
  )
  
  CNa_nose <- 2
  
  # fins
  Xb      <- nose_length + body_length
  Afin    <- 0.5 * (fin_root + fin_tip) * fin_span
  Kfb     <- 1 + (bt_radius / (fin_span + bt_radius))
  
  CNa_fin <- Kfb * (4 * fin_count * (fin_span / bt_diameter)^2) /
    (1 + sqrt(1 + (2 * fin_span / (fin_root + fin_tip))^2))
  
  Xcp_fin <- Xb +
    (fin_sweep / 3) * ((fin_root + 2 * fin_tip) / (fin_root + fin_tip)) +
    (1/6) * (fin_root + fin_tip - (fin_root * fin_tip) / (fin_root + fin_tip))
  
  Cd_fins <- fin_count * 2 * Afin * 0.008 / Aref
  
  # body friction and base drag
  Cd_body <- 0.01
  Cd_base <- 0.029 / sqrt(max(Cd_nose + Cd_body + Cd_fins, 1e-6))
  
  # combined
  Xcp      <- (CNa_nose * Xcp_nose + CNa_fin * Xcp_fin) / (CNa_nose + CNa_fin)
  Cd_total <- Cd_nose + Cd_body + Cd_fins + Cd_base
  
  stability_margin <- (Xcp - cg_measured) / bt_diameter
  
  list(
    Cd               = Cd_total,
    CP               = Xcp,
    stability_margin = stability_margin,
    Cd_nose          = Cd_nose,
    Cd_fins          = Cd_fins,
    Cd_body          = Cd_body,
    Cd_base          = Cd_base
  )
}

# DEFINE FLIGHT SIM

flight_simulation_3d <- function(thrust_curve, prop_mass, dry_mass,
                                 bt_diameter, chute_diameter, chute_delay,
                                 Cd, precision,
                                 wind_speed_ref, wind_dir_deg,
                                 stability_margin,
                                 launch_bearing_deg,
                                 launch_angle_deg) {
  
  if (is.null(thrust_curve) || nrow(thrust_curve) == 0) return(NULL)
  
  thrust     <- approxfun(thrust_curve$time, thrust_curve$thrust,
                          method = "linear", yleft = 0, yright = 0)
  burn_time  <- max(thrust_curve$time)
  m0         <- dry_mass + prop_mass
  bt_area    <- pi * (bt_diameter / 2)^2
  chute_area <- pi * (chute_diameter / 2)^2
  
  total_impulse    <- integrate(thrust,
                                lower = min(thrust_curve$time),
                                upper = max(thrust_curve$time))$value
  exhaust_velocity <- total_impulse / prop_mass
  
  dt <- max(precision, 0.005)
  
  x <- 0; y <- 0; z <- 0
  vx <- 0; vy <- 0; vz <- 0
  angle_rad   <- launch_angle_deg   * pi / 180
  bearing_rad <- launch_bearing_deg * pi / 180
  lx <- sin(angle_rad) * sin(bearing_rad)
  ly <- sin(angle_rad) * cos(bearing_rad)
  lz <- cos(angle_rad)
  m  <- m0
  t  <- 0
  
  t_apogee   <- NA
  chute_open <- FALSE
  wc_gain    <- weathercock_gain(stability_margin)
  
  wind_init    <- wind_velocity(0, wind_speed_ref, wind_dir_deg)
  wx_state     <- wind_init[1]
  wy_state     <- wind_init[2]
  theta        <- 0.5
  k_shape      <- 2.0
  lambda_scale <- wind_speed_ref / gamma(1 + 1/k_shape)
  
  max_steps <- ceiling(1200 / dt)
  out <- data.frame(time     = numeric(max_steps),
                    x        = numeric(max_steps),
                    y        = numeric(max_steps),
                    altitude = numeric(max_steps),
                    velocity = numeric(max_steps),
                    vx       = numeric(max_steps),
                    vy       = numeric(max_steps),
                    vz       = numeric(max_steps))
  i <- 1
  
  repeat {
    
    wind_mean <- wind_velocity(z, wind_speed_ref, wind_dir_deg)
    mu_x      <- wind_mean[1]
    mu_y      <- wind_mean[2]
    mu_mag    <- sqrt(mu_x^2 + mu_y^2)
    
    gust_speed <- rweibull(1, shape = k_shape, scale = lambda_scale)
    wx_state   <- wx_state + theta * (mu_x - wx_state) * dt +
      (gust_speed - mu_mag) * (mu_x / max(mu_mag, 1e-6)) * sqrt(dt)
    wy_state   <- wy_state + theta * (mu_y - wy_state) * dt +
      (gust_speed - mu_mag) * (mu_y / max(mu_mag, 1e-6)) * sqrt(dt)
    wx <- wx_state
    wy <- wy_state
    
    vrel_x   <- vx - wx
    vrel_y   <- vy - wy
    vrel_z   <- vz
    vrel_mag <- max(sqrt(vrel_x^2 + vrel_y^2 + vrel_z^2), 1e-6)
    
    if (chute_open) {
      Cd_eff   <- 1.5
      area_eff <- chute_area
    } else {
      Cd_eff   <- Cd
      area_eff <- bt_area
    }
    
    q      <- 0.5 * rho * vrel_mag^2
    F_drag <- Cd_eff * area_eff * q
    
    Fdx <- -F_drag * (vrel_x / vrel_mag)
    Fdy <- -F_drag * (vrel_y / vrel_mag)
    Fdz <- -F_drag * (vrel_z / vrel_mag)
    
    F_thrust <- thrust(t)
    Fgz      <- -m * g
    
    Fwc_x <- 0; Fwc_y <- 0
    if (is.na(t_apogee) && vz > 0.5) {
      lateral_x <- wx - vx
      lateral_y <- wy - vy
      lat_mag   <- max(sqrt(lateral_x^2 + lateral_y^2), 1e-6)
      coupling  <- wc_gain * q * bt_area * 2.0
      Fwc_x <- coupling * (lateral_x / lat_mag) * min(lat_mag / max(vrel_mag, 1), 1)
      Fwc_y <- coupling * (lateral_y / lat_mag) * min(lat_mag / max(vrel_mag, 1), 1)
    }
    
    ax <- (Fdx + Fwc_x + F_thrust * lx) / m
    ay <- (Fdy + Fwc_y + F_thrust * ly) / m
    az <- (F_thrust * lz + Fdz + Fgz)   / m
    
    if (z <= 0 && t > 0.5 && F_thrust == 0) { ax <- 0; ay <- 0; az <- 0 }
    
    out$time[i]     <- t
    out$x[i]        <- x
    out$y[i]        <- y
    out$altitude[i] <- z
    out$velocity[i] <- sqrt(vx^2 + vy^2 + vz^2)
    out$vx[i]       <- vx
    out$vy[i]       <- vy
    out$vz[i]       <- vz
    i <- i + 1
    
    if (t <= burn_time) {
      m <- max(m - (F_thrust / exhaust_velocity) * dt, dry_mass)
    } else {
      m <- dry_mass
    }
    
    if (!is.na(t_apogee) && !chute_open && (t - t_apogee) >= chute_delay) {
      chute_open <- TRUE
    }
    
    vx <- vx + ax * dt
    vy <- vy + ay * dt
    vz <- vz + az * dt
    x  <- x  + vx * dt
    y  <- y  + vy * dt
    z  <- z  + vz * dt
    t  <- t  + dt
    
    if (is.na(t_apogee) && z > 1 && vz < 0) t_apogee <- t
    if (!is.na(t_apogee) && z <= 0) break
    if (t > 1200 || i >= max_steps) break
  }
  
  out[1:(i-1), ]
}

wind_velocity <- function(altitude_m, wind_speed_ref, wind_dir_deg, z_ref = 10, alpha = 0.14) {
  
  z <- max(altitude_m, z_ref) 
  speed <- wind_speed_ref * (z / z_ref)^alpha
  
  bearing_rad <- (wind_dir_deg * pi / 180)
  Vx <- -speed * sin(bearing_rad) 
  Vy <- -speed * cos(bearing_rad) 
  c(Vx, Vy)
}

# weathercocking
weathercock_gain <- function(stability_margin, CNa_total = 8) {
  gain <- min(stability_margin / 3.0, 1.0) * 0.6
  max(gain, 0)
}

ui <- navbarPage("RRRocket 3D",
                 theme = bs_theme(version = 5, bootswatch = "cerulean"),
                 tabPanel("Setup",
                          page_fluid(
                            navset_card_pill(
                              nav_panel("Rocket", numericInput("dry_mass_kg",  "Dry mass (kg)", value = 0.090),
                                        numericInput("diameter",     "Body tube diameter (mm)", value = 24),
                                        numericInput("body_length", "Body tube length (mm)", value = 300),
                                        numericInput("cg_measured", "CG from nosecone tip (mm)", value = 220)
                              ),
                              nav_panel("Nosecone", selectInput("nose_type", "Nosecone type", choices = c("ogive", "conical", "parabolic")),
                                        numericInput("nose_length", "Nosecone length (mm)", value = 70)
                              ),
                              nav_panel("Fins",
                                        numericInput("fin_count", "Number of fins", value = 3),
                                        numericInput("fin_root", "Root chord (mm)", value = 50),
                                        numericInput("fin_tip", "Tip chord (mm)", value = 25),
                                        numericInput("fin_span", "Semi-span (mm)", value = 30),
                                        numericInput("fin_sweep", "Sweep length (mm)", value = 20)),
                              
                              nav_panel("Chute", 
                                        numericInput("parachute_diameter", "Chute diameter (in)", value = 12)),
                              
                              nav_panel("Launch",
                                        numericInput("rail_length", "Rail length (ft)", value = 3),
                                        numericInput("wind_speed",   "Wind speed (m/s)",     value = 3),
                                        sliderInput("wind_dir", "Where is the wind coming from?", min = 0, max = 360, value = 270),
                                        p("Degrees clockwise from North"),
                                        numericInput("launch_angle",   "Launch angle from vertical (°)", value = 0, min = 0, max = 30),
                                        sliderInput("launch_bearing", "Which direction is the launch angled?",  value = 0, min = 0, max = 360),
                                        p("Also degrees clockwise from North")),
                              
                              nav_panel("Engine",
                                        fileInput("motor_file", NULL, 
                                                  accept = ".eng",
                                                  buttonLabel = "Upload .eng"),
                                        numericInput("parachute_delay",    "Ejection delay (s)", value = 4),
                                        selectInput(
                                          inputId = "engine_choice",
                                          label   = "Or choose an engine:",
                                          choices = c("Select engine..." = "", "A8", "A10", "B4", "B6", "C6", "C11", "D12"),
                                          selected = "B6",
                                          size = 5,
                                          selectize = FALSE
                                        ),
                                        actionButton("visualize_curve", "Visualize Curve")),
                              mainPanel(plotOutput("thrust_curve_plot")))
                          )
                 ),
                 tabPanel("Simulation",
                          sliderInput("precision", "Interval of integration (s)", value = 0.01, min = 0.001, max = 0.1),
                          p("0.01 recommended for accuracy"),
                          actionButton("run", "Simulate", class = "btn-primary")
                 ),
                 tabPanel("Results",
                          plotOutput("altitude_plot"),
                          plotOutput("velocity_plot"),
                          verbatimTextOutput("summary")),
                 
                 tabPanel("3D Track",
                          plotlyOutput("track_3d", width = "600px", height = "600px"),
                          verbatimTextOutput("landing_summary")
                 ),
                 tabPanel("Monte Carlo",
                          sidebarLayout(
                            sidebarPanel(
                              numericInput("mc_runs", "Monte Carlo runs", value = 200, min = 10, max = 1000),
                              sliderInput("montecarlo_precision", "Monte Carlo interval of integration (s)", value = 0.1, min = 0.001, max = 0.5),
                              p("0.1 recommended for time"),
                              sliderInput("chute_delay_std_dev", "Standard deviation of ejection charge delay time (s)", value = 1, min = 0.1, max = 5),
                              p("1 second recommended for Estes motors"),
                              sliderInput("launch_angle_std_dev", "Standard deviation of launch angle (°)", value = 5, min = 0.1, max = 10),
                              p("5 degrees recommended"),
                              sliderInput("wind_speed_std_dev", "Standard deviation of wind speed (%)", value = 5, min = 1, max = 99),
                              p("5 degrees recommended"),
                              sliderInput("wind_dir_std_dev", "Standard deviation of wind direction (°)", value = 30, min = 1, max = 180),
                              p("5 degrees recommended")
                            ),
                            mainPanel(leafletOutput("map", height = 600),
                                      p("Click the map to set launch position"),
                                      actionButton("run_mc", "Run Monte Carlo", class = "btn-warning"),
                                      verbatimTextOutput("landing_pct")
                            )
                          )
                 )
)

parse_thrust_input <- function(motor_file, engine_choice) {
  
  read_eng <- function(lines) {
    lines     <- lines[!grepl("^;", lines)]   # strip comments
    header    <- strsplit(trimws(lines[1]), "\\s+")[[1]]
    prop_mass <- as.numeric(header[5]) / 1000  # grams to kg
    data_lines <- lines[-1]
    pairs     <- lapply(data_lines, function(l) as.numeric(strsplit(trimws(l), "\\s+")[[1]]))
    list(
      thrust_curve = data.frame(time   = sapply(pairs, `[`, 1),
                                thrust = sapply(pairs, `[`, 2)),
      prop_mass    = prop_mass
    )
  }
  
  if (!is.null(motor_file)) {
    lines <- readLines(motor_file$datapath, warn = FALSE)
    read_eng(lines)
  } else {
    req(engine_choice != "")
    path <- file.path("data", paste0(engine_choice, ".eng"))
    req(file.exists(path))
    lines <- readLines(path, warn = FALSE)
    read_eng(lines)
  }
}

server <- function(input, output) {
  
  thrust_curve <- eventReactive(input$visualize_curve, {
    parse_thrust_input(input$motor_file, input$engine_choice)$thrust_curve
  })
  
  results <- eventReactive(input$run, {
    parsed <- parse_thrust_input(input$motor_file, input$engine_choice)
    aero   <- compute_aero(
      nose_type   = input$nose_type,
      nose_length = input$nose_length / 1000,
      body_length = input$body_length / 1000,
      bt_diameter = input$diameter    / 1000,
      fin_count   = input$fin_count,
      fin_root    = input$fin_root    / 1000,
      fin_tip     = input$fin_tip     / 1000,
      fin_span    = input$fin_span    / 1000,
      fin_sweep   = input$fin_sweep   / 1000,
      cg_measured = input$cg_measured / 1000
    )
    sim <- flight_simulation_3d(
      thrust_curve       = parsed$thrust_curve,
      dry_mass           = input$dry_mass_kg,
      prop_mass          = parsed$prop_mass,
      bt_diameter        = input$diameter / 1000,
      chute_diameter     = input$parachute_diameter * 0.0254,
      chute_delay        = input$parachute_delay,
      Cd                 = aero$Cd,
      precision          = input$precision,
      wind_speed_ref     = input$wind_speed,
      wind_dir_deg       = input$wind_dir,
      stability_margin   = aero$stability_margin,
      launch_angle_deg   = input$launch_angle,
      launch_bearing_deg = input$launch_bearing
    )
    list(sim = sim, aero = aero)
  }) 
  
  landing_pct_val <- reactiveVal(NULL)
  
  monte_carlo <- eventReactive(input$run_mc, {
    parsed <- parse_thrust_input(input$motor_file, input$engine_choice)
    aero   <- compute_aero(
      nose_type   = input$nose_type,
      nose_length = input$nose_length / 1000,
      body_length = input$body_length / 1000,
      bt_diameter = input$diameter    / 1000,
      fin_count   = input$fin_count,
      fin_root    = input$fin_root    / 1000,
      fin_tip     = input$fin_tip     / 1000,
      fin_span    = input$fin_span    / 1000,
      fin_sweep   = input$fin_sweep   / 1000,
      cg_measured = input$cg_measured / 1000
    )
    n <- input$mc_runs
    landings <- lapply(1:n, function(i) {
      sim <- flight_simulation_3d(
        thrust_curve = parsed$thrust_curve,
        prop_mass = parsed$prop_mass,
        dry_mass = input$dry_mass_kg,
        bt_diameter = input$diameter / 1000,
        chute_diameter = input$parachute_diameter * 0.0254, #in to m conversion
        chute_delay = max(0,   rnorm(1, input$parachute_delay, input$chute_delay_std_dev)),
        Cd = aero$Cd,
        precision = input$montecarlo_precision,
        wind_speed_ref = max(0, rnorm(1, input$wind_speed, 0.01 * input$wind_speed * input$wind_speed_std_dev)),
        wind_dir_deg = rnorm(1, input$wind_dir, input$wind_dir_std_dev),
        stability_margin = aero$stability_margin,
        launch_angle_deg = max(0,   rnorm(1, input$launch_angle, input$launch_angle_std_dev)),
        launch_bearing_deg = input$launch_bearing
      )
      tail(sim, 1)[, c("x", "y")]
    })
    do.call(rbind, landings)
  })
  
  output$thrust_curve_plot <- renderPlot({
    req(thrust_curve())
    tc <- thrust_curve()
    
    thrust_fn     <- approxfun(tc$time, tc$thrust, yleft = 0, yright = 0)
    total_impulse <- integrate(thrust_fn, min(tc$time), max(tc$time))$value
    max_thrust <- max(tc$thrust)
    
    ggplot(tc, aes(time, thrust)) +
      geom_area(fill = "darkseagreen4", alpha = 0.68) +
      geom_line() +
      geom_hline(yintercept = 0, color = "gray") +
      annotate("text",
               x = max(tc$time), y = max(tc$thrust),
               label = sprintf("Total impulse: %.2f Ns", total_impulse),
               hjust = 1, vjust = 1, size = 6, fontface = "bold", color = "black") +
      annotate("text",
               x = max(tc$time), y = max(tc$thrust) - (max(tc$thrust) / 10),
               label = sprintf("Max thrust: %.2f N", max_thrust),
               hjust = 1, vjust = 1, size = 6, fontface = "bold", color = "black") +
      scale_y_continuous(limits = c(0, max(tc$thrust + (max(tc$thrust)/7)))) +
      labs(x = "time (s)", y = "thrust (N)", title = "Engine thrust vs. time") +
      theme_classic() +
      theme(plot.title = element_text(face = "bold", size = 16))
  })
  
  output$altitude_plot <- renderPlot({
    req(results())
    ggplot(results()$sim, aes(time, altitude * m_to_ft)) +
      geom_line() +
      geom_hline(yintercept = 0, color = "gray") +
      labs(x = "time (s)", y = "altitude (ft)", title = "Altitude") +
      theme_classic() +
      theme(plot.title = element_text(face = "bold", size = 16))
  })
  
  output$velocity_plot <- renderPlot({
    req(results())
    r <- results()$sim
    req(nrow(r) > 0)
    ggplot(r, aes(time, vz * m_to_ft)) +
      geom_line() +
      geom_hline(yintercept = 0, color = "gray") +
      labs(x = "time (s)", y = "z velocity (ft/s)", title = "Vertical velocity") +
      theme_classic() +
      theme(plot.title = element_text(face = "bold", size = 16))
  })
  
  output$track_3d <- renderPlotly({
    req(results())
    r <- results()$sim
    
    plot_ly(r, x = ~x * m_to_ft, y = ~y * m_to_ft, z = ~altitude * m_to_ft,
            type = "scatter3d", mode = "lines",
            line = list(color = ~altitude * m_to_ft, colorscale = c('#86CEFA', '#003396'), width = 5)) |>
      add_trace(x = tail(r$x, 1) * m_to_ft,
                y = tail(r$y, 1) * m_to_ft,
                z = 0,
                type = "scatter3d", mode = "markers",
                marker = list(color = "red", size = 6),
                name = "landing") |>
      layout(scene = list(
        xaxis = list(title = "East (ft)"),
        yaxis = list(title = "North (ft)"),
        zaxis = list(title = "Altitude (ft)")
      )
      )
  })
  
  output$summary <- renderPrint({
    req(results())
    r    <- results()$sim
    aero <- results()$aero   
    cat("apogee ", round(max(r$altitude) * m_to_ft), "ft\n")
    cat("time at apogee ", round(r$time[which.max(r$altitude)], 3), "s\n")
    cat("max velocity ", round(max(r$velocity) * m_to_ft, 3), "ft/s\n")
    cat("total flight time ", round(max(r$time), 3), "s\n")
    cat("stability margin ", round(aero$stability_margin, 2), "calibers\n")
  })
  

  output$landing_pct <- renderPrint({
    req(monte_carlo())
    lc    <- monte_carlo()
    drift <- sqrt(lc$x^2 + lc$y^2) * m_to_ft
    
    cat(sprintf("95th percentile distance from launchpad: %.0f ft\n", quantile(drift, 0.95)))
    cat(sprintf("Max distance from launchpad: %.0f ft\n", max(drift)))
    
    if (!is.null(landing_pct_val())) {
      cat(sprintf("%.1f%% of simulated rockets land in the polygon", landing_pct_val()))
    }
  })
  
  output$mc_plot <- renderPlotly({
    req(monte_carlo())
    lc <- monte_carlo()
    
    plot_ly(lc, x = ~x * m_to_ft, y = ~y * m_to_ft,
            type = "scatter", mode = "markers",          # scatter not lines
            marker = list(color = "steelblue", size = 6, opacity = 0.5)) |>
      layout(
        xaxis = list(title = "East drift (ft)", range = c(0, max(x))),
        yaxis = list(title = "North drift (ft)", range = c(0, max(y))),
        title = "Monte Carlo landing scatter"
      )
  })
  
  launch_point <- reactiveVal(list(lat = 38.89, lng = -77.03))
  
  observeEvent(input$map_click, {
    launch_point(list(lat = input$map_click$lat, lng = input$map_click$lng))
  })
  
  output$map <- renderLeaflet({
    leaflet() |>
      addProviderTiles("Esri.WorldImagery") |>
      setView(lat = 36.9052, lng = -81.0768, zoom = 5) |>
      addDrawToolbar(
        polylineOptions     = FALSE,
        circleOptions       = drawCircleOptions(),
        markerOptions       = FALSE,
        circleMarkerOptions = FALSE,
        rectangleOptions    = drawRectangleOptions(),
        polygonOptions      = drawPolygonOptions(),
        editOptions         = FALSE
      )
  })
  
  observe({
    lp <- launch_point()
    leafletProxy("map") |>
      clearGroup("launch") |>
      addMarkers(lng = lp$lng, lat = lp$lat,
                 label = "Launch pad",
                 group = "launch")
  })
  
  drawn_polygon <- reactiveVal(NULL)
  
  observeEvent(input$map_draw_new_feature, {
    drawn_polygon(input$map_draw_new_feature)
  })
  
  observeEvent(monte_carlo(), {
    req(drawn_polygon())
    req(launch_point())          # make sure launch point is set
    
    lc     <- monte_carlo()
    lp     <- launch_point()     # <-- use the reactiveVal, not input$launch_marker
    lat0   <- lp$lat
    lng0   <- lp$lng
    
    lc$lat <- lat0 + (lc$y / 111320)
    lc$lng <- lng0 + (lc$x / (111320 * cos(lat0 * pi / 180)))
    
    # safely extract polygon coordinates
    poly_coords <- drawn_polygon()$geometry$coordinates[[1]]
    poly_mat    <- do.call(rbind, lapply(poly_coords, function(p) c(p[[1]], p[[2]])))
    
    # close the polygon ring if not already closed
    if (!identical(poly_mat[1,], poly_mat[nrow(poly_mat),])) {
      poly_mat <- rbind(poly_mat, poly_mat[1,])
    }
    
    polygon_sf  <- st_polygon(list(poly_mat)) |> st_sfc(crs = 4326)
    points_sf   <- st_as_sf(lc, coords = c("lng", "lat"), crs = 4326)
    inside      <- st_within(points_sf, polygon_sf, sparse = FALSE)[,1]
    pct         <- mean(inside) * 100
    landing_pct_val(pct)
    
    lc$color <- ifelse(inside, "green", "red")
    
    leafletProxy("map") |>
      clearMarkers() |>
      clearGroup("launch") |>        # re-add launch marker after clearMarkers wipes it
      addCircleMarkers(lng = lc$lng, lat = lc$lat,
                       radius = 4, color = lc$color,
                       fillOpacity = 0.6, stroke = FALSE) |>
      addMarkers(lng = lng0, lat = lat0,
                 label = "Launch pad",
                 group = "launch") |>
      addPopups(lng = lng0, lat = lat0,
                popup = sprintf("%.1f%% land inside safe zone", pct))
  })
  
}

shinyApp(ui = ui, server = server)