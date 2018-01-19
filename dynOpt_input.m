function param = dynOpt_input(varargin)
% the input parser is inside here

    % create and define the parser structure
    % if the matlab function is compiled, can only accept strings as param
    % due to command line limitations
    p = inputParser;
    p.FunctionName = 'param parser for dynOpt';
    p.CaseSensitive = true;
    
        % REQUIRED
        addRequired(p,'elecKE');
            % (average) kinetic energy of the electron, eV
        % specify a range in the form of a$b for the GA/parameter sweep
        addRequired(p,'space_matrix');
            % [total repeat units in system] n*[type {other propreties}]
            % little rectangles: -1 + [n rectangles], a.u.
            % sine waves: -2 + n*[coeff + phase], coeff in a.u.
        addRequired(p,'time_matrix');
            % [type {other properties}]
            % piecewise: -1 + [0 rise 1 fall phase], a.u.
                % phase from 1-4 (off rise on fall)
            % sine waves: -2 + tau + n*[coeff + phase], a.u.
        addRequired(p,'Lx');
            % length in nm
        addRequired(p,'Nx');
            % number of grid points - should be power of 2 
            
        % OBSERVABLES SWITCHES
        addParameter(p,'doNorm',0,@(x) any(x==[0 1]));
            % 1: <Psi|Psi> 
            % 2: <Psi|Psi> in forbidden regions (doEnergy must be 1)
        addParameter(p,'doMV',0,@(x) any(x==[0 1]));
            % 1: <v>
        addParameter(p,'doMV2',0,@(x) any(x==[0 1]));
            % 1: <v^2>
        addParameter(p,'doX',0,@(x) any(x==[0 1]));
            % 1: <x>
        addParameter(p,'doX2',0,@(x) any(x==[0 1]));
            % 1: <x^2>
        addParameter(p,'doFlux',0,@(x) any(x==[0 1]));
            % 1: calculate the flux out of the flux points
        addParameter(p,'doEnergy',0,@(x) any(x==[0 1 2 3]));
            % 0: off
            % 1: check T only
            % 2: check T + V
        % Combination observables switches
        addParameter(p,'doEff',0,@(x) any(x==[0 1]));
            % 1: return efficiency as the fitness value to be optimized
        addParameter(p,'doDiffusion',0,@(x) any(x==[0 1]));
            % 0: off
            % 1: calculate diffusion by 2nd moments (rev. mod. phys 1983)
            
        % SIMULATION
        addParameter(p,'doAvg',0,@(x) any(x==[0 1]));
            % 0: dynFitness returns full output
            % 1: dynFitness returns one value (averaged)
        addParameter(p,'doGA',0,@(x) any(x==[0 1]));
            % 0: don't do it
            % 1: genetic algorithm search
        addParameter(p,'doSawtooth',0,@(x) any(x==[0 1]));
            % 0: unrestricted optimization
            % 1: all harmonics restricted to 1/2 max amplitude of first
        addParameter(p,'doSweep',0);
            % 0: don't do it
            % >1: number of points in the linspace of the sweep
        addParameter(p,'doMovie',0,@(x) x>0);
            % from end going backwards, visualize [x]% of the simulation
            % [% 1:full(rho) 2:diag(rho)]
        addParameter(p,'doHybrid',0,@(x) any(x==[0 1 2]));
            % 0: off
            % 1: fmin
            % 2: patternsearch (recommended since gradients not easily estimated)
        addParameter(p,'saveIntermediate',0,@(x) any(x==[0 1]));
            % save the result in a file everytime dynFitness is called
        addParameter(p,'verbose',0,@(x) any(x==[0 1]));
            % print timing information to stdout
        addParameter(p,'ode_solver',113,@(x) any(x==[23 45 113]));
            % choose different nonstiff solvers (113 is the best)
        addParameter(p,'continue',0,@(x) x>0);
            % continue time for propagation
            
        % AB INITIO parameters
        addParameter(p,'prop',1,@(x) any(x==[1 2 3 4]))
            % 1: split operator
            % 2: short iterative lanczos
            % 3: chebyshev
            % 4: density operator
        addParameter(p,'cutoff',0,@(x) 0<=x && x<=1)
            % cut off first [x]% of the time series before observables
        addParameter(p,'time_multiplier',20,@(x) x >~ 1)
            % multiple of largest Vt time constant to simulate for
        addParameter(p,'x0',0);
            % center of the wavepacket, nm
        addParameter(p,'sigmax',10); 
            % 3-sigma spread of WP, nm
        addParameter(p,'bias',5);
                % add a (-) to get a thermal ensemble at the temperature
                % specified by elecKE (rule of thumb: 25 meV ~ 300 K) -
                % this really only has meaning when paired with density
                % matrix propagation (prop = 4)
            % gaussian
                % 1 = running to the left only
                % 2 = spreading left and right
                % 3 = right bias
            % plane wave
                % 4 = left moving plane wave
                % 5 = left + right
                % 6 = right moving
            % random phases
                % 7 = random phase uniform amplitude WF
                % 8 = gaussian envelope, random phase
            % eigenfunctions
                % 9.fraction = eigenfunction, f is the fraction of the 
                % potential
        addParameter(p,'doRatchet',1,@(x) any(x==[0 1]));
            % 1: turn on time dependent potential
        addParameter(p,'doNIP',0,@(x) any(x==[0 1]));
            % 0: turn off the absorbing potential (a.k.a. periodic system)
            % 1: turn on the absorbing potential
        addParameter(p,'doBath',0,@(x) any(x(1)==[0 1 2]));
            % [type {bath parameters}]
            % 0: No bath, prop 1-4 valid
            % 1: (S)SH, prop 2-3 valid [(some parameters here)]
            % 2: Lindblad, only prop 4 valid [T gamma]
            % 3: Lindblad, pure-decoherence [T gamma]
        addParameter(p,'effMass',1,@(x) x>0);
            
    % parse the params and insert them into the param structure
    parse(p,varargin{:});
    
    % take params out of p
    param = p.Results;
    
    % DECLARE PHYSICAL CONSTANTS (ATOMIC UNITS)
        % from http://physics.nist.gov/cuu/index.html
        param.Eh = 27.21138505; % hartree energy in eV
        param.a0 = 5.2917721092E-11; % bohr radius (m), to convert nm 
        param.auE = 5.14220652E11; % a.u. of electric field, V/m
        param.k = 8.6173324E-5; % eV/K
        param.aut = 2.418884326502E-17; % a.u. of time
    
    % set up the job title
    param.job_title = strtok(getenv('PBS_JOBID'),'.');
    param.tor = datestr(now,'dd-mm-yy_HH-MM-SS');

    % flux point - % of the grid point number
    param.fPoint = floor(0.90*param.Nx);

    % Set the first element to 999 if you want to use an elliptical
    % 'transmissionless' absorber (set = 0 to turn off absorber)
    % 777 left (bottom) side onLy
    % 888 right (top) side onLy
    param.absorbx = [999 floor(0.95*param.Nx) param.Nx];
    
    if mod(log2(param.Nx),1) ~= 0
        warning('Nx is not a power of 2. FFTW will not be at its most efficient!')
        fftw('planner', 'patient');
    end
    
    % set up doGA or doSweep
    if param.doGA
        % for GA, each fitness evaluation returns efficiency
        param.doAvg = 1;
        param.doMV = 1;
        param.doMV2 = 1;
        param.doEff = 1;
    elseif param.doSweep
        % for a sweep, each run must have <v>, <v^2>, <KE>
        param.doMV = 1;
        param.doMV2 = 1;
        param.doEnergy = 2;
    end
    
    %% bath parameters

    if param.doBath(1) == 1
        % note that bit-ordering begins from 0 (MATLAB starts at 1), as well as
        % numbering of the spinor components. So bath modes run from 0 -> Nbath-1,
        % and spinor components run from 0 -> 2^Nbath - 1

        % spectral density (polynomial or some function)
        % J = eta * omega * e^(omega/omega_c)
        param.Je = @(x) 1*x;

        % number of bath modes
        param.bathN = 3;

        % sample the spectral density at discrete energy steps (linear for now)
        param.e_sample = linspace(0.5,1.3,param.bathN+1)/param.Eh;

        % density of states
        param.Je_dos = zeros(1,param.bathN);
        for a = 1:param.bathN
            param.Je_dos(a) = 1/(param.e_sample(a+1) - param.e_sample(a));
        end

        % discrete spectral density
        param.Je_disc = param.Je(param.e_sample(1:end-1));

        % bath interaction term
        param.Je_coupling = sqrt(param.Je_disc./param.Je_dos);

        % number of spinor components (# of wavefunctions)
        spinorComps = 2^param.bathN;

        % build the relaxation operators
        % diagonal energy matrix for the number operator (initialized as vector)
        epsilon = zeros(spinorComps,1);
        for a = 1:spinorComps
            for b = 1:param.bathN
                if bitget(a-1,b) == 1 % if the b'th bit (b'th mode) of the spinor component is set
                    epsilon(a) = epsilon(a) + param.e_sample(b); % add the energy of the mode to that spinor component
                end
            end
        end

        % find the correct spinor components that when bitxor with the Je_coupling
        % (0:bathN-1) index yields the current spinor component looked at
        relax_index = zeros(spinorComps,param.bathN);
        %Dkl = zeros(spinorComps);
        for a = 1:spinorComps % spinor being looked at
            for b = 0:param.bathN-1
                for c = 1:spinorComps
                    if bitxor(2^b,c-1) == a-1 % i = bitxor(2^k,j), psi_i = couple_k * psi_j
                        relax_index(a,b+1) = c;
                        %Dkl(a,c) = param.Je_coupling(b+1);
                    end
                end
            end
        end

        % build the dephasing operator
        % global dephasing parameter
        param.cbar = 1;

        % inelastic bias - smaller -> smaller coefficients
        param.sigma_epsilon = 0.01;

        % build the reduced dephasing index
        if mod(param.bathN,2) == 0
            % maximum size of the dephasing matrices
            dephase_size = (param.bathN/2)^2;
        else
            dephase_size = (param.bathN-1)/2 + ((param.bathN-1)/2)^2;
        end
        dephase_index = zeros(spinorComps,dephase_size);
        cij = zeros(spinorComps,dephase_size);
        %Okl = zeros(spinorComps);
        pairs = nchoosek(0:param.bathN-1,2);
        for a = 1:spinorComps % spinor being looked at
            i_ = 1;
            for b = 1:spinorComps % look over all other spinor comps
                if sum(bitget(a-1,1:param.bathN)) == sum(bitget(b-1,1:param.bathN)) % same # of excitations
                    for c = 1:nchoosek(param.bathN,2) % look over possible cij combinations
                        if bitxor(a-1,b-1) == sum(2.^pairs(c,:)) % bitxor(i,j) == sum(2.^[k l])
                            dephase_index(a,i_) = b;
                            cij(a,i_) = param.cbar*exp(-(param.e_sample(pairs(c,1)+1) - param.e_sample(pairs(c,2)+1)).^2/(2*param.sigma_epsilon^2))/(param.bathN*(param.bathN-1));
                            %Okl(a,b) = param.cbar*exp(-(param.e_sample(pairs(c,1)+1) - param.e_sample(pairs(c,2)+1)).^2/(2*param.sigma_epsilon^2))/(param.bathN*(param.bathN-1));
                            i_ = i_ + 1;
                        end
                    end
                end
            end
        end

        param.epsilon = reshape(epsilon,1,[]); param.relax_index = relax_index;
        param.cij = cij; param.dephase_index = dephase_index;
        clear epsilon relax_index a b c i_ dephase_index cij dephase_size
    end
    
    if param.doBath(1) >= 2
        % temperature (kT - in eV)
        param.kT = param.k*param.doBath(2)/param.Eh;
        
        % friction
        param.gamma = param.doBath(3);
    end
end
