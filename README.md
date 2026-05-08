# RRRocket 3D
Model rocket flight simulator built with R Shiny.
https://tatecommission.shinyapps.io/rrrocket3d/

# Run locally
install.packages(c("shiny","bslib","plotly","ggplot2","leaflet","leaflet.extras"))
shiny::runApp()

# Why RRRocket?
Easier to use than competitors with intuitive UI and easily customizable paramaters.
You don't need to know how to code, and you don't need to know complex aerodynamics.

One-of-a-kind Monte Carlo launch zone selection tool lets you assess launch site viability
combining satellite imagery with Monte Carlo simulated trajectories. For example, what is 
the simulated probability of a rocket landing in your given area with user-defined uncertainty
paramaters?

Convienent web format. Some of the complex aerodynamics are sacrificed in order for the app
to run smoothly as an R Shiny app on the internet.

## Physics
Models kinematics by integrating the kinematic equations with Euler's Method using 
user-defined time interval. Motor mass decreases as motor burns.
Uses Barrowman equations with Galejs body-lift correction to calculate center of pressure,
which changes as the motor spends fuel and the rocket changes position.
Wind is distributed according Weibull distribution, modified by altitude according to power law.
Both constant winds and sudden gusts are accounted for with this technique.
Atmosphere modeled with full 4-layer ISA standard up to 80km as a foundation for future updates.

## Current limitations
Not yet accurate at higher altitudes
Looking to implement more accurate stability simulation with 6 degrees of freedom, rivaling OpenRocket
No staging or dual parachutes yet
Euler integration can accumulate error, considering implementing diff eq