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
<sup>+</sup>: These authors contributed equally to this manuscript.

This repository serves as the home for the open-source Peanut Butter and Joints biomimetic robotic grasper. The design, modelings, and characterization of this grasper is detailed in the manuscript references above. This repository contains the design files (SolidWorks part and assembly files, STLs for 3D printed components, and engineering drawings (limited) ), simulation files and associated code, and experimental data associated with the above manuscript. 

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

### Directory: [Arduino Controller]()
This directory contains a folder entitled "Code Files" which hold all files required to implement the position controller described in the manuscript. 
 - "DCMotor.cpp" and "DCMotor.h" are library files which provide code to interface with the motor that actuates the hand.
 - "handrobot.ino" is the main Arduino/Teensy code file that implements the controller and can be used to conduct the ball drop tests.
 
### Directory: [CAD]()
This directory contains the SolidWorks part files, STL files, and engineering drawings associated with the grasper and its various subcomponents. Each of these folders is further organized by the subsystem to which the files belong.
 - "Hand:" files associated with the 3D printed hand, including the palm, phalanges, and encoder mounts.
 - "Damper:" files associated with the concentric ring viscous damper (peanut butter not included)
 - "Hand Stand:" files associated with the experimental test stand used to folder the grasper in place during drop tests.
 - "Ball Stand:" files associated with the experimental test stand use to hold the ball servo and laser distance sensor used during ball drop tests.
 - "Misc:" miscellanious files associated with other components of the hand (e.g., molds for casting silicone components) and experimental setups (e.g., the pendulum joints).
Note: not all sub-components appear in all three subdirectories.

### Directory: [Modeling and Experimental Data]()
This directory contains all files and code associated with the modeling and characterization of the grasper and its subcomponents.
