function dynOpt_func(varargin)
% minimize the function dynFitness
% M.U.T.H.E.R. - dispatches mcFitness or dynFitness

    % for compiled function param, need to convert from string to numbers
    for a = 1:nargin
        temp = [];
        if ischar(varargin{a})
            % str2num only converts if input is valid matlab syntax for a
            % number or matrix, leaving the name of name-value pairs as
            % strings
            temp = str2num(varargin{a});
        end
        if ~isempty(temp)
            varargin{a} = temp;
        end
    end

    %% check the first argument of varargin and parse it appropriately
    
    % lower and upper bounds for the GA and sweeper
    % for AbI, the first block holds the space variables to be optimized, 
    % while the second block is the time vars
    % for MC, the boundaries basically go between 0 and 1, since the
    % actual scattering gradient shape is normalized to be between 0
    % and 1

    % linear inequality constraints
    A = [];
    bineq = [];
    i_ = 1;
    
    if strcmp(varargin{1},'AbI')
        fprintf('\n*****************************\n');
        fprintf('AB INITIO')
        fprintf('\n*****************************\n')
        param = dynOpt_input(varargin{2:end});
        f = fullfile('defectOutput',['dynOpt_' param.job_title '_' param.tor]);
        save(f,'param');
        
        [param.space_matrix,lb{1},ub{1}] = split_input(param.space_matrix);
        [param.time_matrix,lb{2},ub{2}] = split_input(param.time_matrix);
        lb = cell2mat(reshape(lb,1,[]));
        ub = cell2mat(reshape(ub,1,[]));
        
        % number of variables that we are looking at, by GA or sweep
        nvars = sum(isnan(param.space_matrix)) + sum(isnan(param.time_matrix));
        
        % number of -1 or -2 (tile types)
        param.unique_tiles = sum(param.space_matrix<0);

        % location of tile designators - add 1 to the end marker since to take
        % slices out of the matrix we go from ind -> next ind - 1
        param.space_ind = [find(param.space_matrix < 0) length(param.space_matrix)+1];
        param.time_ind  = [find(param.time_matrix < 0)  length(param.time_matrix)+1];
        
        % work up the space matrix
        if param.doSawtooth
            for a = 1:param.unique_tiles
                % ms = matrix slice, cut out each unique tile from the space_matrix
                % to examine it
                ms = param.space_matrix(param.space_ind(a):param.space_ind(a+1)-1);
                switch ms(1)
                    case -1
                        % rectangles
                        num_nan = sum(isnan(ms));
                        for b = 1:num_nan-1
                            % insert pairs of [-1 1], each step must be
                            % smaller than the previous one for a sawtooth
                            temp = zeros(1,nvars);
                            temp(i_) = -1;
                            temp(i_+1) = 1;
                            A = [A ; temp];
                            i_ = i_ + 1;
                        end
                        i_ = i_ + 1; % jump index since the for loop runs from i_ -> i_+1
                    case -2
                        for b = 2:length(ms) % n*[coeff phase]
                            if isnan(ms(b))
                                % if nan is given (i.e. optimize)
                                if mod(b,2) == 0
                                    % coefficient
                                    if b == 2 % first harmonic
                                        a_1 = i_;
                                        i_ = i_ + 1;
                                    else % other harmonics max 1/2
                                        temp = zeros(1,nvars);
                                        temp(a_1) = 1/(b/2);
                                        temp(i_) = 1;
                                        A = [A ; temp];
                                        i_ = i_ + 1;
                                    end
                                else
                                    % phase factor
                                    % sawtooth does not affect phase, skip
                                    i_ = i_ + 1;
                                end
                            end
                        end
                end
            end
        end
        
        % create the anonymous function
        FitnessFunction = @(x) dynFitness(x,param);
    elseif strcmp(varargin{1},'MC')
        fprintf('\n*****************************\n');
        fprintf('MONTE CARLO')
        fprintf('\n*****************************\n')
        
        param = mcInput(varargin{2:end});
        % save a copy of the parameters
        f = fullfile('defectOutput',['dynOpt_' param.job_title '_' param.tor]);
        save(f,'param');
        
        % parse ranges into NaNs and generate upper/lower boundaries
        [param,lb,ub] = mcVary(param,[]);
        
        % number of variables
        nvars = length(ub);
        
        % create the anonymous function
        FitnessFunction = @(x) mcFitness(x,param);
    else
        error('You did not specify something dynOpt can do!')
    end
    
    % set the b 
    if param.doSawtooth
        bineq = zeros(size(A,1),1);
    end
    
    %% set up and run the GA    
    if param.doGA
        fprintf('\n*****************************\n');
        fprintf('GENETIC ALGORITHM')
        fprintf('\n*****************************\n')

        % Start with the default options
        options = gaoptimset;

        % Modify options setting
        options = gaoptimset(options,'PopInitRange', [lb;ub]);
        options = gaoptimset(options,'Display','iter');
        options = gaoptimset(options,'UseParallel','always');
        options = gaoptimset(options,'PopulationSize',2*str2double(getenv('PBS_NP'))); % twice number of processors

        switch param.doHybrid
            case 0
                % no hybrid functional
            case 1
                % fmin
                hybridopts = optimoptions('fmincon','Display','iter','UseParallel','always');
                options = gaoptimset(options,'HybridFcn',{@fminunc,hybridopts}); 
            case 2
                % pattern search
                hybridopts = psoptimset;
                hybridopts = psoptimset(hybridopts,'Display','iter');
                hybridopts = psoptimset(hybridopts,'CompletePoll','on'); % poll all points on mesh
                hybridopts = psoptimset(hybridopts,'UseParallel','always');
                options = gaoptimset(options,'HybridFcn',{@patternsearch,hybridopts});
        end

        % shuffle the rng seed
        rng('shuffle','simdTwister')

        % parallel pool initialization
        parallelOpen

        % run the genetic algorithm
        [x,fval,exitflag,ga_output,population,score] = ga(FitnessFunction,nvars,A,bineq,[],[],lb,ub,[],[],options);

        % parallel pool shutdown
        parallelClose

        % calculate <<v>> and <<KE>>
        param.doMV2 = 0;
        param.doEnergy = 1;
        param.doEff = 0;
        output = dynFitness(x,param);
        meanv = output(1);
        meane = output(2);

        save(f,'meane','meanv','x','fval','exitflag','ga_output','population','score','options','-append');
    end
    
    %% do a parameter sweep    
    if param.doSweep
        
        fprintf('\n*****************************\n');
        fprintf('PARAMETER SWEEP')
        fprintf('\n*****************************\n')
        
        assert(nvars==length(param.doSweep),...
            'Not enough sweep ranges were not given for $ variables!')
        
        % NaN's indicate parameters to do a sweep, create an array that
        % goes from left to right, each dimension is a NaN that is swept
        % the lb and ub that is generated at the beginning of this
        % function will be reused here
        if isscalar(param.doSweep)
            sweep_data = cell(param.doSweep,1);
        else
            sweep_data = cell(param.doSweep);
        end
        doSweep = param.doSweep;
        vary_ls = cell(numel(doSweep),1);
        for a = 1:numel(vary_ls)
            vary_ls{a} = linspace(lb(a),ub(a),doSweep(a));
        end

        % parallel pool initialization
        parallelOpen
        
        parfor a = 1:prod(doSweep)
            % convert the linear index to subscripts
            vary_ind = cell(numel(doSweep),1);
            [vary_ind{:}] = ind2sub(doSweep,a);
            vars = cellfun(@(x,y) x(y),vary_ls,vary_ind);
            
            % actually do the calculation
            sweep_data{a} = FitnessFunction(vars);
        end
        
        % parallel pool shutdown
        parallelClose
        
        save(f,'sweep_data','lb','ub','-append');
    end
end