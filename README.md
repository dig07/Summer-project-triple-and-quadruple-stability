# Summer-project-triple-and-quadruple-stability
Student: Diganta Bandopadhyay<br>
Supervisor: Dr Silvia Toonen<br><br>
Summer project by undergraduate student investigating the dynamical stability of hierarchical triple and quadruple stellar systems. Stellar systems are 
picked from probability distributions of orbital parameters; orbital seperations, eccentricities and mass ratios. Main purpose of study is to expoose any
underlying stability biases within orbital parameter distributions. <br><br>
<b>Triple stellar systems key results</b>
<ul>
  <li>Orbits with eccentricities picked from a thermal distribution
  are more likely to undergo interactions between the inner
  binary and the triple star.</li>
  <li>T14 orbital seperation distribution produces more systems with binary-triple
  interactions than OBin.</li>
  <li>Flat and 1/q mass distributions both produce similar numbers of stable stellar systems.</li>
</ul>
<b>Quadruple stellar systems key results</b>
<ul>
  <li>Massive stars are observed to have more companions and tend
  to be in more compact orbits (observational evidence). There
  may be some link between the number of companions and the
  size of orbit, in the quadruple case the inner binary orbital separation (stable) is
  pushed to lower orbital separations (more often) than in the
  triple case.</li>
</ul>
<br>

<b>Simulation Process</b><br>
For each simulation: 
<ul>
 <li>Initial_PDFs.py file creates a .txt file containing all the stellar systems to be checked for dynamical stability. Disbtributions for orbital 
 parameters are specified in this file.</li>
 <li> Simulate.py imports the .txt file created by Initial_PDFs.py containing all the systems to be checked. All systems are checked for stability. 
 Outputs text file containing counter arrays, these are arrays which keep track of the number of stable systems with each of the orbital elements picked 
 in the Initial.py file.</li>
 <li>Analysis.py imports in the text files created by Initial_PDFs.py and Simulate.py and plots the parent and stable histograms across orbital parameters.</li>
</ul><br>

The Simulate.py file has to check stability across ~10^12 stellar systems to make sure there is enough data points to produce histograms with clear trends. 
To do this within reasonable execution time, we code the intensive part of the process in Fortran95 this is the simulation.f95 file which is compiled into a 
python library (simulation.so) using f2py. Additionally the multiprocessing library is used to evaluate the stability of systems in parallel, the simulation script
is most efficent if run on a computer with multiple CPUs, (the number of processes should not exceed the number of logical CPU cores). 
