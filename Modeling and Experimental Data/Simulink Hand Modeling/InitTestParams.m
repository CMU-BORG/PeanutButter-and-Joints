g=9.81;
Axle2AxleDist_mm = 30;
Axle2AxleDist_TIP_mm = 35;
JointDampingScaling=1;
MotorTorque = 0.4;

S1.tau_s = [-0.071,0.145,-0.154,0.029]; %polynomial coefficients.  note that this is based on the angle of the PP
S1.K = [1.02 -0.54 0.45]; % polynomial coefficients. Based on the difference between the current angle and the previous angle
S1.JointDamping = 0.0142*JointDampingScaling;
S1.InitAngle=0;
S1.lowerLimit=-20;
S1.upperLimit=90;
S1.Axle2AxleDist_mm = 30;
S1.ContactStiffness = 1e6;
S1.ContactDamping=1e3;

S2.tau_s = [0.056 0.016 -0.132 0.015]; %polynomial coefficients.  note that this is based on the angle of the the MP relative to the PP
S2.K = [1.06 -0.76 0.40]; % polynomial coefficients. Based on the difference between the current angle and the previous angle
S2.JointDamping = 0.0105*JointDampingScaling;
S2.InitAngle=0;
S2.lowerLimit=-20;
S2.upperLimit=90;
S2.Axle2AxleDist_mm = 30;
S2.ContactStiffness = 1e6;
S2.ContactDamping=1e3;

S3.tau_s = [-0.103 0.102 -0.052 -0.019]; %polynomial coefficients.  note that this is based on the angle of the DP relative to the MP
S3.K = [0.38 -0.09 0.13]; % polynomial coefficients. Based on the difference between the current angle and the previous angle
S3.JointDamping = 0.0081*JointDampingScaling;
S3.InitAngle=0;
S3.lowerLimit=-20;
S3.upperLimit=90;
S3.Axle2AxleDist_mm = 30;
S3.ContactStiffness = 1e6;
S3.ContactDamping=1e3;

Palm.ContactStiffness = 1e6;
Palm.ContactDamping=1;

Ball.X_Pos = 0;
Ball.Y_Pos = 300; %height above palm in mm
Ball.Z_Pos = 0;
Ball.R=33; %millimeters
Ball.Mass=56; %grams

CollisionOpacity=0.1;

