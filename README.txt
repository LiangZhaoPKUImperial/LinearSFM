===============================================================================================
Linear SFM: A Hierarchical Approach to Solving Structure from Motion Problems by Decoupling the Linear and Nonlinear Components
Version: 1.0
===============================================================================================

Copyright (C) 2015 Liang Zhao, Shoudong Huang, and Gamini Dissanayake
University of Technology, Sydney, Australia

C/C++ sourse code for Linear SFM: A Hierarchical Approach to Solving Structure from Motion Problems by Decoupling the Linear and Nonlinear Components

Authors:  Liang Zhao         -- liang.zhao@imperial.ac.uk 
          Shoudong Huang     -- Shoudong.Huang@uts.edu.au
	  Gamini Dissanayake -- Gamini.Dissanayake@uts.edu.au

	  Hamlyn Centre for Robotic Surgery
          Department of Computing
          Faculty of Engineering
          Imperial College London, United Kingdom          

	  Centre for Autonomous Systems
          Faculty of Engineering and Information Technology
          University of Technology, Sydney
          NSW 2007, Australia
-----------------------------------------------------------------------------------------------
License
-----------------------------------------------------------------------------------------------

Linear SFM: by Liang Zhao, Shoudong Huang, Gamini Dissanayake is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

1) You can freely use and modify this code.

2) If you want to distribute code based on this one, it has to be done under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

If you use this code for academic work, please reference:

      Liang Zhao, Shoudong Huang, and Gamini Dissanayake,
      Linear MonoSLAM: A Linear Approach to Large-Scale Monocular SLAM Problems,
      IEEE International Conference on Robotics and Automation (ICRA), 2013.     

3) For commercial uses, please contact the authors.

-----------------------------------------------------------------------------------------------
Quick start
-----------------------------------------------------------------------------------------------
Please run the code on x64 platform.

On Windows platform, directly open the "LinearSFM.sln" solution and uncomment some codes corresponding to one dataset in _tmain() function.

On Linux platform, install and run ParallaxBA as fowllows:

- sudo apt-get install cmake
- sudo apt-get install libeigen3-dev
- sudo apt-get install libsuitesparse-dev

- cd linux
- mkdir build
- cd build
- cmake ..
- make
- sudo make install

- cd DataForC
- LinearSFM -path RS468_C -p Pose_RS468.txt -f Feature_RS468.txt -num 466 -type Monocular
- LinearSFM -path NC3500_C -p Pose_NC3500.txt -f Feature_NC3500.txt -num 3499 -type Stereo
- LinearSFM -path RS90_C -p Pose_RS90.txt -f Feature_RS90.txt -num 88 -type Monocular

-----------------------------------------------------------------------------------------------
Support
-----------------------------------------------------------------------------------------------

Questions, comments, suggestions, discussions and bug reports are welcomed. 

Please email to liang.zhao@imperial.ac.uk

