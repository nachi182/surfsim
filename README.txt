      IBM Research, 2014
      Martin Suchara
      SURFACE CODE SIMULATOR
----------------------------------
----------------------------------
SHORT DESCRIPTION:
Surface code threshold simulator that assumes that syndrome measurements can fail. Parameters allow simulation of both the phenomenological and circuit error model. Also simulates leakage. This is an optimized version that doesn't match on the complete graph but rather constructs a graph where for each vertex we find the nearest neighbor in the "up", "down", "right", and "left" direction and only consider edges inside of this rectangle.
The first line specifies the lattice sizes that the simulator evaluates. The second line specifies the number of repetitions of the simulation for each lattice size. For example, 100000 means that a random error is generated 100000 times and the simulator attempts to correct the error. The output is the percentage of the generated errors that were corrected successfully. The third line swithces between the circuit and phenomenological error model.
The following three lines specify how is the depolarizing error in each iteration generated. For example, pMin = 0.00, pMax = 0.01, pStep = 0.0005 means that a series of simulations is performed where the error is generated using the depolarizing error model starting with p=0% and following to p=1.0% by increasing the error in 0.05% increments.
The following two lines specify the two-qubit errors in the circuit model and the measurement error, respectively. For example, p2 = 10 means that each of the 15 two-qubit errors has probability 10 * p/15, and q = 8 means that the syndrome measurement fails with probability 8 * p.
The last line determines if a leakage reduction is attempted.
The next four lines encode the leakage parameters that are simulated. pLeakUp is the probability of leaving the computational state. pLeakDown is the probability of a leaked qubit returning back to the computational state. pleakUp = 0.5 means that the probability of leakage of a gate is 0.5 * p and pLeakDown = 0.1 means that the leakage reduction probability of a gate is 0.1 * p.
The simulation parameters are specified in in/parameters.txt. The file has the following format:
The simulation result is saved in the out directory in separate files for each leakage parameter and lattice size. A new directory is created with a name that encodes the leakage paramaters that are simulated. If the directory already exists it is overwritten. The out directory also contains scripts that visualize the data. For visualization, copy the .txt files to the appropriate directories and type ./get_plot to produce a .ps and .pdf file with the visualizations.
To compile type make clean; make. To run the simulation type ./qecsim.
XSize = 5 7 9
circuitErr = yes
correctLeak = no
iterations = 1000
p2 = 1
pLeakDown = 1.0
pLeakUpMax = 1.0
pLeakUpMin = 1.0
pLeakUpStep = 0.1
pMax = 0.01
pMin = 0.00
pStep = 0.0005
q = 1
1.1 COMPILATION
1.2 SIMULATION PARAMETERS
1.3 SIMULATION OUTPUT
