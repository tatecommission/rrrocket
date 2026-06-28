# RRRocket 3D

RRRocket 3D is a model rocket flight simulator built in R Shiny. The goal is a tool that is both physically honest and easy enough to use that you do not need an engineering background to get something meaningful out of it.

## What makes it different

Most simulators give you an apogee and a flight profile. RRRocket 3D adds a Monte Carlo landing predictor on top of that. After you set up your rocket and run a nominal simulation, you can run hundreds of perturbed flights in the background and see where they all land on a real satellite map. You can then draw a polygon over your field and get an instant readout of what percentage of simulated flights land inside it. That makes it genuinely useful for pre-launch safety planning, not just performance estimation.

The 3D trajectory viewer also colors the flight by phase: red during motor burn, yellow through coast, and pale green during parachute descent. This makes it immediately obvious how much the wind is pushing the rocket during each part of the flight, which a flat altitude-versus-time plot does not show you.

## Physics

Most of the aerodynamics are modeled off of and sometimes simplified from "The Practical Calculation of the Aerodynamic Characteristics of Slender Finned Vehicles" by James Barrowman and the OpenRocket technical documentation by Sampo Niskanen. The rocket is split into three main parts: nosecone, body tube, and fins. The drag coefficients for each are calculated separately in alignment with Niskanen.

Drag is decomposed into four components. Nose pressure drag, which is the drag caused by the nosecone pushing air out of the way during flight, uses the half-angle sine-squared formula and depends primarily on nosecone type and geometry. Skin friction drag, the drag caused by air rubbing against all exterior surfaces, is computed from the Reynolds-number-based turbulent flat-plate formula (Barrowman eq. 3.78), switching to a roughness-limited value (eq. 3.80) above the critical Reynolds number, with a Mach compressibility correction applied both subsonic and supersonic (eqs. 3.82-3.83). Base drag, caused by the low-pressure zone that forms directly below the rocket, uses the Hoerner Mach-dependent formula: 0.12 + 0.13M² for M < 1, 0.25/M for M ≥ 1, applied at each timestep. Parasite drag refers to all drag besides base drag, and all components of parasite drag are assumed to scale roughly equally with Mach number.

The atmosphere follows the ISA troposphere model (288.15 K at sea level, -6.5 K/km lapse rate). Wind uses a 1/7-power-law altitude shear profile with an Ornstein-Uhlenbeck turbulence process, and the weathercocking effect is accounted for. Flight is integrated with Euler's method on a user-defined time interval. As the fuel burns, the center of gravity and stability margin update at each timestep.

## Monte Carlo

The Monte Carlo mode perturbs ejection delay, launch angle, wind speed and direction, dry mass, drag coefficient, and propellant mass using random samples from Gaussian distributions with user-defined standard deviations. Each run returns only the landing coordinate to keep things fast even at 200 or more runs. Results are plotted on a satellite map, and users can draw a polygon over their intended landing area to compute what fraction of simulated flights land inside it.

## Engine data

Motors use the standard RASP .eng format. Built-in curves cover common Estes A through G motors. For anything larger, the .eng file can be downloaded from [thrustcurve.org](https://www.thrustcurve.org) and uploaded directly.