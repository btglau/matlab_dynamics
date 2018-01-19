# matlab_dynamics
A repository of MATLAB code that performs quantum dynamics on a grid, as well as semiclassical Monte-Carlo for semiconductors. The purpose of this code is to investigate the time-dependent velocity / position as a function of various potentials. The quantum dynamics is done on a grid, while the semiclassical Monte-Carlo uses a continuous axis. The quantum dynamics module is more for toy potentials, while the Monte-Carlo module is capable of including realistic scattering mechanisms found in semicondcutors.

Quantum dynamics uses fft/ifft (i.e. Fourier grid Hamiltonian), and can do:
- Time-dependent potentials, arbitrary shapes through piecewise / Fourier series
- Conservative: split-operator, short-iterative lanczos, chebyshev
- Dissipative: stochastic surrogate hamiltonian, density matrix + lindblad
https://journals.aps.org/pre/abstract/10.1103/PhysRevE.93.062128

Semiclassical Monte-Carlo (Newtonian analytic propagation + QM scattering rates):
- Phonon, photon, and impurity scattering
http://onlinelibrary.wiley.com/doi/10.1002/aenm.201701000/full

Both modules output results in a MATLAB structure, and have an argument parser for accepting input. The main function which calls dynamics / aggregates output, dynOpt_func, can do full parameter space sweeps or use MATLAB's built in stochastic search techniques (such as genetic algorithms) to search parameter space.

All code is "heavily" MATLAB optimized, but does not take advantage of the newest R2016 feature that rolls bsxfun functionality into element-by-element operations (i.e. bsxfun(@times,A,B) is now A.*B).

### Required MATLAB toolboxes:
MATLAB

Optimization Toolbox

Image Processing Toolbox

Statistics and Machine Learning Toolbox

Global Optimization Toolbox

Parallel Computing Toolbox

MATLAB Distributed Computing Server
