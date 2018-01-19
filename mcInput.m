function param = mcInput(varargin)
% the input parser is inside here

    % create and define the parser structure
    % if the matlab function is compiled, can only accept strings as param
    % due to command line limitations
    p = inputParser;
    p.FunctionName = 'param parser for manchego carlo';
    p.CaseSensitive = true;
    
        % RANGE PARAMETERS - a r b
        % the GA will optimize within that range or a parameter sweep will 
        % be done within that range
        addParameter(p,'comp_x',0); % compositional gradient
        addParameter(p,'comp_y',0); % compositional gradient
        addParameter(p,'comp_z',0); % compositional gradient
        addParameter(p,'imp_x',0); % impurity gradient
        addParameter(p,'imp_y',0); % impurity gradient
        addParameter(p,'imp_z',0); % impurity gradient
            % rectangles: -1 + [repeat units] + n*[0 -> 1]
            % sine waves: -2 + [repeat units] + n*[coeff + phase]
            % polynomial: -3 + [repeat units] + [p1*x^n .. pn*x^1 + pn+1]
            % sawtooth: -4 + [repeat units] + [%]
            % nothing: 0
        addParameter(p,'Ex',[0 0 0]);
        addParameter(p,'Ey',[0 0 0]);
        addParameter(p,'Ez',[0 0 0]);
            % [V/cm V/cm Hz]: DC-offset AC-amplitude Hz
        addParameter(p,'kick',[0 0 0]);
            % [mu(eV) sigma(eV) rate(Hz)]
            % a bath that adds energy to electrons at a certain rate, given
            % by an average and spread (normally distributed)
        addParameter(p,'photon',[0 0]);
            % [eV J/m^3]
            % photon bath (infrared light) - supercedes the kicking
            % mechanism if both are specified
        addParameter(p,'comp_range',[0 0.5]);
        addParameter(p,'imp_range',[0 1]);
            % if a fourier series / linear gradient is used, its range
        addParameter(p,'Lx',1000);
        addParameter(p,'Ly',1000);
        addParameter(p,'Lz',1000);
            % length PER PERIOD in nm, default 1 micron
        addParameter(p,'T0',300);
            % crystal temperature, kelvin
        addParameter(p,'Te',300);
            % (average) kinetic energy of the e-, eV
            % used as the 'kT' variable in P(E) = f(E)g(E)
        addParameter(p,'ni',1E16);
            % highest doping density in the system (in cm^-3)
        addParameter(p,'co',0)
            % single value of strain (between 0 and 1, in %)
        addParameter(p,'specularity',0.5)
            % 1: purely specular reflection
            % 0: purely diffusive reflection
            
        % OBSERVABLES SWITCHES
        addParameter(p,'doX',0,@(x) any(x==[0 1]));
            % 1: <x>
        addParameter(p,'doX2',0,@(x) any(x==[0 1]));
            % 1: <x^2>
        addParameter(p,'doMV',0,@(x) any(x==[0 1]));
            % 1: <v>
        addParameter(p,'doMV2',0,@(x) any(x==[0 1]));
            % 1: <v^2>
        addParameter(p,'doFlux',0,@(x) any(x==[0 1]));
            % 1: calculate the flux out of the flux points
        addParameter(p,'doEnergy',0,@(x) any(x==[0 1 2 3]));
            % 1: KE
        addParameter(p,'doValley',0,@(x) any(x==[0 1]));
            % 1: calculate average occupation of lower valley
        addParameter(p,'doMovie',0,@(x) any(x==[0 1]));
            % 0->1: movie of x, y, and z coordinates starting @ time
            % collects sim_slice number of frames
            
        % Combination observables switches
        addParameter(p,'doEff',0,@(x) any(x==[0 1]));
            % 1: return efficiency as the fitness value to be optimized
        addParameter(p,'doDiffusion',0,@(x) any(x==[0 1]));
            % 0: off
            % 1: calculate diffusion by 2nd moments (rev. mod. phys 1983)
            
        % SIMULATION
        addParameter(p,'doAvg',0,@(x) any(x==[0 1]));
            % 0: mcFitness returns full output
            % 1: mcFitness returns one value (averaged)
        addParameter(p,'doGA',0,@(x) any(x==[0 1]));
            % 0: don't do it
            % 1: genetic algorithm search
        addParameter(p,'doSawtooth',0,@(x) any(x==[0 1]));
            % 0: unrestricted optimization
            % 1: all harmonics restricted to 1/2 max amplitude of first
        addParameter(p,'doSweep',0);
            % 0: don't do it
            % >1: number of points in the linspace of the sweep (per $)
        addParameter(p,'doHybrid',0,@(x) any(x==[0 1 2]));
            % 0: off
            % 1: fmin
            % 2: patternsearch (recommended since gradients not easily estimated)
        addParameter(p,'saveIntermediate',0,@(x) any(x==[0 1]));
            % save the result in a file everytime dynFitness is called (for GA)
        addParameter(p,'verbose',0,@(x) any(x==[0 1]));
            % print timing information to stdout
        addParameter(p,'sim_time',1E-9,@(x) x > 0)
            % simulation time, s, default 1 ns
        addParameter(p,'sim_slice',1E3,@(x) x > 1)
            % number of times during sim to collect observables
        addParameter(p,'dt_slice',1,@(x) x > 1)
            % slice the free flights into smaller free flights (for
            % calcuations with electric fields)
        addParameter(p,'doPBCx',0,@(x) x>0);
        addParameter(p,'doPBCy',0,@(x) x>0);
        addParameter(p,'doPBCz',0,@(x) x>0);
            % periodic boundary conditions (corresponding E should be 0)
        addParameter(p,'Ne',5000,@(x) x > 1)
            % number of electrons in simulation
        
    % parse the params and insert them into the param structure
    parse(p,varargin{:});
    
    % take params out of p
    param = p.Results;
    
    % DECLARE PHYSICAL CONSTANTS (ATOMIC UNITS)
        % from http://physics.nist.gov/cuu/index.html
        param.Eh = 27.21138505; % eV/Ha
        param.a0 = 5.2917721092E-11; % m/a0
        param.auE = 5.14220652E9; % a.u. of electric field, V/cm
        param.k = 8.6173324E-5; % eV/K
        param.aut = 2.418884326502E-17; % s/a.u. of time
        param.erm = 9.10938356E-31; % kg/me, electron rest mass
        param.wEv = 6.2415091E18; % eV/J
        param.f = 8.23872336E-8; % Eh/a0, N
    
    % set up the job title
    param.job_title = strtok(getenv('PBS_JOBID'),'.');
    param.tor = datestr(now,'dd-mm-yy_HH-MM-SS');
    
    % make sure some conditions are met
    % set up doGA or doSweep
    if param.doGA
        % for GA, each fitness evaluation returns efficiency
        param.doAvg = 1;
        param.doMV = 1;
        param.doMV2 = 1;
        param.doEff = 1;
    elseif param.doSweep
        % <v>, <v^2>, <KE>
        param.doMV = 1;
        param.doMV2 = 1;
        param.doEnergy = 1;
    end
    % need the values of q to get average values of <q^2>
    if param.doMV2
        param.doMV = 1;
    end
    if param.doX2
        param.doX = 1;
    end
    % don't do movie if averaging is on
    if param.doAvg
        param.doMovie = 0;
    end
end