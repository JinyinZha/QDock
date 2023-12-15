  Pose sampling in molecular docking is burdensome for classical computers. Quantum computers (QC) are possible hardware-solution to computational expensive problems. However, QCs are usually critical to input. Therefore, pose sampling in molecular docking could not directly be solved by QCs.<p>
  QDock is a python package encoding pose sampling in molecular docking for quantum conputers. In brief, it translates pose sampling into QC-intepretable QUBO models. There are two methods in QDock package, including "Grid Point Matching" and "Feature Atom Matching". Grid Point Matching is close to the sampling performance of Glide SP, while Feature Atom Matching is computationally more efficient. A brief illustration of both methos are given below. QDock also integrates a QC-simulator to solve the sampling problem.<p><p>
Recommended Configuration<p>
CentOS Platform<br>
Python   3.9<br>
PyQUBO   1.2.0<br>
Prody    2.2.0<br>
Numpy    1.23.4<br>
Scipy    1.7.1<br>
<br>
AutoDockFR should be downloaded and installed. Its "bin" path should be added to enviromental variables. AutoDockFR is at https://ccsb.scripps.edu/adfr/downloads/ <p>
If you found QDock useful, please cite the followiung article.<p>
https://doi.org/10.1021/acs.jctc.3c00943
