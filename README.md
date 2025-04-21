# PeanutButter-and-Joints
Paper: "PB&J: Peanut Butter and Joints for Damped Articulation"

An open-source, low-cost, biomimetic robotic grasper platform modeled on the human hand, incorporating the effects of viscoelastic joint mechanics.

An archived version of this repository is available through Zenodo: [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.XXXXXXXXX.svg)]()

Authors:

* Avery S. Williamson<sup>+</sup>
* Michael J. Bennington<sup>+</sup>
* Ravesh Sukhnandan<sup>+</sup>
* Mrinali Nakhre
* Yeumin Mao
* Victoria A. Webster-Wood

<i>All authors are affiliated with Carnegie Mellon University</i>

+: These authors contributed equally to this manuscript.

This repository serves as the home for the open-source Peanut Butter and Joints biomimetic robotic grasper. The design, modeling, and characterization of this grasper is detailed in the manuscript references above. This repository contains the design files (SolidWorks part and assembly files, STLs for 3D printed components, and engineering drawings (limited) ), simulation files and associated code, and experimental data associated with the above manuscript. 

> [!NOTE]
> If you are interesting in utilizing or modifying this grasper or any of its subcomponents for your own <i>in roboto</i> investigations of human like grasping or robotic grasper development, please cite using the following DOI:
> [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.XXXXXXXXX.svg)]() 
> <b>Above needs to be updated with correct DOI for either Living Machines paper or preprint archive</b>
>
> Williamson, Bennington, and Sukhnandan, <i>et al.</i>, "PB&J: Peanut Butter and Joints for Damped Articulation." (2025). <i>arXiv</i> <b>Place doi here</b>


## Overview:

The files in the repository are organized into three subdirectories:

- Arduino Controller

- CAD

- Modeling and Experimental Data

### Directory: [Arduino Controller](https://github.com/CMU-BORG/PeanutButter-and-Joints/tree/main/Arduino%20Controller)
This directory contains a folder entitled "Code Files" which hold all files required to implement the position controller described in the manuscript. 
 - "DCMotor.cpp" and "DCMotor.h" are library files which provide code to interface with the motor that actuates the hand.
 - "handrobot.ino" is the main Arduino/Teensy code file that implements the controller and can be used to conduct the ball drop tests.
 
### Directory: [CAD](https://github.com/CMU-BORG/PeanutButter-and-Joints/tree/main/CAD)
This directory contains the SolidWorks part files, STL files, and engineering drawings associated with the grasper and its various subcomponents. Each of these folders is further organized by the subsystem to which the files belong.
 - "Hand:" files associated with the 3D printed hand, including the palm, phalanges, and encoder mounts.
 - "Damper:" files associated with the concentric ring viscous damper (peanut butter not included)
 - "Hand Stand:" files associated with the experimental test stand used to folder the grasper in place during drop tests.
 - "Ball Stand:" files associated with the experimental test stand use to hold the ball servo and laser distance sensor used during ball drop tests.
 - "Misc:" miscellanious files associated with other components of the hand (e.g., molds for casting silicone components) and experimental setups (e.g., the pendulum joints).
Note: not all sub-components appear in all three subdirectories.

### Directory: [Modeling and Experimental Data](https://github.com/CMU-BORG/PeanutButter-and-Joints/tree/main/Modeling%20and%20Experimental%20Data)
This directory contains all files and code associated with the modeling and characterization of the grasper and its subcomponents. This directory is broken down by the associated model/characterization experiments performed. There are subfolders for:
 - Concurrent Flexion testing
 - Couette Flow Damper Modeling and the associated design parameter sweeps
 - Pendulum Damper Modeling and the associated pendulum drop test videos and Tracker files
 - Simulink Hand Modeling 
 
#### Sub-directory: [Concurrent Flexion](https://github.com/CMU-BORG/PeanutButter-and-Joints/tree/main/Modeling%20and%20Experimental%20Data/Concurrent%20Flexion) 
File List:
 - CompareFlexion.m: MATLAB script to perform range-of-motion normalized comparisons between finger flexion with and without a parallel elastic element. This script also calculates the associated pair-wise correlation coefficients reported in the manuscript.
 - Tracking Data Flexion.xlsx: position tracking data from videos of flexion experiments with and without a parallel elastic element.
 
#### Sub-directory: [Couette Flow Damper Modeling](https://github.com/CMU-BORG/PeanutButter-and-Joints/tree/main/Modeling%20and%20Experimental%20Data/Couette%20Flow%20Damper%20Modeling) 
File List:
 - Effective_Damping_Coeff_v2.m: MATLAB function that calculated the effective damping coefficient (G) for a given input geometry for the concentric cylinder damper
 - Parameter_Sweeps_v2.m: MATLAB script that performs the parametric sweeps of the damper design space and calculates the heatmaps of effective damping coefficients.
 - PlotFins.m: MATLAB function to plot the cross-section of the concentric ring dampers for visualization.
 
#### Sub-directory: [Pendulum Damper Modeling](https://github.com/CMU-BORG/PeanutButter-and-Joints/tree/main/Modeling%20and%20Experimental%20Data/Pendulum%20Damper%20Modeling)
Folder List:
 - 20250414_Parameter: folder containing the final distributions for the friction and damping parameters resulting from the bootstrapping method.
 - No Damper Pendulum Tests: folder containing the videos, Tracker files, and resulting positional data from the pendulum drop tests without a viscous damper
 - Peanut Butter Damper Pendulum Tests: folder containing the video of all drop experiments, the Tracker file, and resulting data. All experiments are within the single video
 
File List:
 - AllTrialError.m: MATLAB function that calculates the error between the pendulum model and a set of experimental data for a given parameter set. Used to calculate the objective function in the optimization.
 - dynamics.m: MATLAB function that calculated the model angular acceleration for a given position, velocity, and set of parameters.
 - PB_J__Parameter_Estimation.pdf: Text document explaining the details of the parameter estimation routine.
 - Pendulum_ParameterEstimation.m: MATLAB script that conducts the bootstrap parameter estimation of the frictional and damping parameters from the pendulum drop tests. To rerun the Monte Carlo simulations to obtain the average mdoel results, set the "if_save_friction" and "if_save_damping" parameters to 0. This will reload the saved parameters from the folder listed above. If you want to rerun the optimization, set these variables to 1. <b>Important!:</b> Before you run the optimization, create a new save folder and update the "outfolder" variables or else you will overwrite the results. Additionally, the optimization is slow and you should expect it to take many hours. It is best to let it run overnight.
 - theta_model.m: MATLAB function that returns a vector of model angles for a given set of parameters and initial conditions. These values are calculated by integrating the presented dynamical model and linearly reinterpolating at a set of time points that align with the experimental data time points.
 
#### Sub-directory: [Simulink Hand Modeling](https://github.com/CMU-BORG/PeanutButter-and-Joints/tree/main/Modeling%20and%20Experimental%20Data/Simulink%20Hand%20Modeling)
Folder List:
 - CAD: holds the CAD files of the hand used as the geometry of the hand model
 - slprj: autogenerated folder from Simulink.
 
File List: 
 - CatchingBall.avi/.mp4: video files showing the simulation successfully capturing a ball. The kinematics of these simulations is shown in the manuscript.
 - CheckFingerFunction.slx: Simulink model utilizing Simscape Multibody to simulate the dynamics of the hand, the forces transmission of the pulleys, and the contact interactions with the ball.
 - InitTestParams.m: MATLAB script that initializes the simulation paramters for the Simulink model. This script must be run before the Simulink model can be run.
 - PlotResults.m: Script to plot the results of the Simulink simulations and generate the figure presented in the manuscript. This script can be run after the Simulink simulation is complete.