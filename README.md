# TractionController
ProjectV4 contains ME675 Course Project, Nonlinear MPC Traction controller. Some notes:
  - This project was very difficult and did not get successful results until a week before it was due
  - All that to say, this repo is a mess.
  - Good luck
      - Go to ECE675 Project/ProjectV4/
      - Model_v5.m assembles the model using an instance of Car.m. CarWithDefaults.m inherits and sets default values. Model_v5 references many other scripts which serve other functions, and programmatically builds the differential equations for the vehicle model. This is built with MATLAB symbolic variables. If you want to simulate a different vehicle, adjust any of the parameters of the car.m instance.
      - .m simulates the vehicle under split traction conditions. The vehicle starts on dry pavement (mu=0.9 on all 4 tires), and eventually hits ice on left hand tires (mu=0.2) with high traction maintained on right hand tires.
      - The 5th order DAE is converter to CasADi symbolic variables, so that the CasADi implementation of the ADAS solver can be used for iterative optimization. This is done by BuildIntegrator.m, and MPC formulation is done by BuildMPC.m
      - This read me is getting too long, message me if you want more detail. its called ECE675 but the course was actually ME675. oopsy typo
