# 2D-vertex-dynamics-model

****** From Make to Run ********  
rm CMakeCache.txt  
cmake .  
make  
./2dv  
*****************************  


***************************  
Ignore the following comments after simulation runs  

Line tension by PCP = 0  
Fluctuation = 0  
Power of inner product (PCP) = 2  
Pulse Period = 55  
Phase Shift along x-axis = 0.7854  
Phase Shift along y-axis = 1.5708  
Random seed = 20110412  
****************************  

****************************  
_parameters.h  
PERIOD_PARAVIEW: vtk file is output per PERIOD_PARAVIEW (100000 prefered)  
TIME_CELL_DIVISION: cell cycle duration  
****************************  

****************************  
2dv.cpp  
at LINEs 560, 561  
 axis.x = cos(theta);  
 axis.y = sin(theta);  
Cell division axis is given by these two lines  
****************************  

 
