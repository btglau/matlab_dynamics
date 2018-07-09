# matlab_dynamics

*Note*: Since it will be quite hard to figure out how to use this code, the hope is that the propagate / mcFitness functions, which contain the code for the dynamics, can be used as skeletons for other work.

A repository of MATLAB code that performs quantum dynamics on a grid, as well as semiclassical Monte-Carlo for semiconductors. The purpose of this code is to investigate the time-dependent velocity / position as a function of various potentials. The quantum dynamics is done on a grid, while the semiclassical Monte-Carlo uses a continuous axis. The quantum dynamics module is more for toy potentials, while the Monte-Carlo module is capable of including realistic scattering mechanisms found in semicondcutors.

Quantum dynamics uses fft/ifft (i.e. Fourier grid Hamiltonian), and can do:
- Time-dependent potentials, arbitrary shapes through piecewise / Fourier series
- Conservative: split-operator, short-iterative lanczos, chebyshev
- Dissipative: stochastic surrogate hamiltonian, density matrix + lindblad
- Paper: https://journals.aps.org/pre/abstract/10.1103/PhysRevE.93.062128
- https://doi.org/10.1016/0021-9991(91)90137-A

I would like to highlight propagate_1D_lvn.m in particular - it makes extensive use of nested functions with ode45/113 to simplify passing things around. I do not know if this is best practice, or even best efficiency ... please let me know!

Semiclassical Monte-Carlo (Newtonian analytic propagation + QM scattering rates):
- Phonon, photon, and impurity scattering
- Paper: http://onlinelibrary.wiley.com/doi/10.1002/aenm.201701000/full
- https://journals.aps.org/rmp/pdf/10.1103/RevModPhys.55.645

Both modules output results in a MATLAB structure, and have an argument parser for accepting input. The main function which calls dynamics / aggregates output, dynOpt_func, can do full parameter space sweeps or use MATLAB's built in stochastic search techniques (such as genetic algorithms) to search parameter space.

All code is "heavily" MATLAB optimized, but does not take advantage of the newest R2016 feature that rolls bsxfun functionality into element-by-element operations (i.e. bsxfun(@times,A,B) is now A.*B).

### Required MATLAB toolboxes:
- MATLAB
- Optimization Toolbox
- Image Processing Toolbox
- Statistics and Machine Learning Toolbox
- Global Optimization Toolbox
- Parallel Computing Toolbox
- MATLAB Distributed Computing Server

### Getting started
Check example input.

Run from within MATLAB with:

dynFitness [commands] \
or \
dynFitness(Position,Name,Value,...)

or, via command line (https://www.mathworks.com/help/matlab/ref/matlablinux.html) \
matlab -r [commands] \
or \
matlab -nodisplay -nosplash -r [commands]

Input file -> dynOpt_func -> selects quantum dynamics or monte carlo -> processes input with dynOpt_input or mcInput -> calls appropriate dynamics function dynFitness / mcFitness; dynFitness further calls the 'initialize' family functions to set up the grid, then calls various 'propagate' functions.
