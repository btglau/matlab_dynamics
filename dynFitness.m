function [ fit_value ] = dynFitness(vars,param)
    %% tic
    totalStart = tic;
    
    fprintf('\n==========================\n')
    fprintf('DYNAMIC FITNESS!!!\n')
    fprintf('==========================\n')

    %% assign the GA variables
    % vars is what the GA sends in - from 1 : number of NaN's in space matrix
    % is the potential variables, and from that + 1 to the end is time
    var_counter = sum(isnan(param.space_matrix));
    param.space_matrix(isnan(param.space_matrix)) = vars(1:var_counter);
    param.time_matrix(isnan(param.time_matrix)) = vars(var_counter+1:end);
    
    % param to atomic units
    param.elecKE = param.elecKE/param.Eh;
    param.Lx = 1E-9*param.Lx/param.a0;
    param.x0 = 1E-9*param.x0/param.a0;
    param.sigmax = 1E-9*param.sigmax/param.a0; 
        % 3x this value is 3sigma (nm)
        
    % define grid spacing (i.e. where the grid basis f'n are)
    % grid is from 0 <= x <= L, with sampling L/N
    % in the fourier method, the function is sampled y_n = y(nL/N) with n =
    % 0 -> N-1
    param.dx = param.Lx/param.Nx;
    
    % define the grid
    param.X = (0:param.Nx-1).'*param.dx;

    % momentum grid point spacing pi/gridLength determines the highest K
    % allowed.
    param.dkx = (2*pi)/param.Lx;
    
    %% Initialize the operators and set up some values for the wavefunction prior to running dynamics

    % Initialize the kinetic matrix
    [param.T,param.kxs,param.Tr] = initializeT_1D(param);
    if param.verbose
        fprintf('Finished kinetic matrix initialization: ')
        toc(totalStart)
    end

    % Initialize the ratchet potential
    [param.V,param.imagV] = initializeV_1D(param);
    if param.verbose
        fprintf('Finished potential intialization: ')
        toc(totalStart)
    end
    
    % Initialize time parameters
    param = initializeTime_1D(param);
    if param.verbose
        fprintf('Finished time parameter intialization: ')
        toc(totalStart)
    end
    
    % initialize WF
    Psi = initializeWF_1D(param,param.elecKE);
    if param.verbose
        fprintf('Finished wavefunction initialization: ')
        toc(totalStart)
    end

    %% Propagation and observable collection

    if param.verbose
        fprintf('Normalization before prop %f',sum(abs(Psi(:)).^2))
    end
        if param.prop == 4
            [Psi,output] = propagate_1D_lvn(Psi,param);
        else
            [Psi,output] = propagate_1D(Psi,param);
        end        
    if param.verbose
        if param.prop == 4 % density matrix
            fprintf('\nNormalization after prop: %f\n',trace(Psi))
        else
            fprintf('\nNormalization after prop: %f\n',sum(abs(Psi(:)).^2))
        end
        fprintf('Finished propagation: ')
        toc(totalStart)
    end        
    
    %% output the GA values
    
    if param.doAvg
        % if doAvg is enabled, any observable is returned as an average
        
        % the actual simulation time (1-cutoff)
        simul_time = param.simul_time * (1-param.cutoff);
        tau_max = max(param.all_taus);
        
        % average over integer multiples of a period (if ratcheting)
        if param.doRatchet
            % find how many periods (tau_max) fit in the total simulation time
            per_length = floor(simul_time/tau_max) * tau_max;
            if param.prop == 4
                per_end = find(output.time_vector >= param.simul_time*param.cutoff + per_length,1,'first');
                simul_time = output.time_vector(per_end) - output.time_vector(1);
            else
                per_end = round(per_length/param.dt);
            end
        else
            if param.prop == 4
                % set the averaging to the end of the time_vector
                per_end = length(output.time_vector);
            else
                % yes, Nt*(1-cutoff) is proper, but this emulates the
                % behaviour in propagate_1D to be numerically safe
                per_end = param.Nt - round(param.Nt*param.cutoff);
            end
        end
        
        if param.doMV
            % <p>/m
            if param.prop == 4
                output.MV = trapz(output.time_vector(1:per_end),output.MV(1:per_end))/simul_time;
            else
                output.MV = mean(output.MV(1:per_end));
            end
        end

        if param.doMV2     
            % <p^2>/m^2
            if param.prop == 4
                output.MV2 = trapz(output.time_vector(1:per_end),output.MV2(1:per_end))/simul_time;
            else
                output.MV2 = mean(output.MV2(1:per_end));
            end
        end

        if param.doX
            % <x>
            if param.prop == 4
                output.MX = trapz(output.time_vector(1:per_end),output.MX(1:per_end))/simul_time;
            else
                output.MX = mean(output.MX(1:per_end));
            end
        end

        if param.doX2
            % <x^2>
            if param.prop == 4
                output.MX2 = trapz(output.time_vector(1:per_end),output.MX2(1:per_end))/simul_time;
            else
                output.MX2 = mean(output.MX2(1:per_end));
            end
        end

        if param.doFlux
            % fluxxx
            % not averaging, just a straight up time integral
            if param.prop == 4
                output.flux_right = trapz(output.time_vector(1:per_end),output.flux_right(1:per_end));
                output.flux_left = trapz(output.time_vector(1:per_end),output.flux_left(1:per_end));
            else
                output.flux_right = param.dt*trapz(output.flux_right(1:per_end));
                output.flux_left = param.dt*trapz(output.flux_left(1:per_end));
            end
        end

        if param.doEnergy > 0
            % average energy (0, 1, or 2 columns)
            if param.prop == 4
                output.en = trapz(output.time_vector(1:per_end),output.en(1:per_end,:),1)/simul_time;
            else
                output.en = mean(output.en(1:per_end,:),1);
            end
        end
        
        % fit_value is what is returned to the GA; the temporary value
        % ('temp') accumulates other values from the propagation and can 
        % spit it into an file that tracks observales at intermediate times
        fit_value = [output.MV output.MV2 output.MX output.MX2 output.flux_right output.flux_left output.en];
        
        if param.doEff
            % efficiency defined as <v>^2/<v^2>
            fit_value = output.MV^2/output.MV2;

            % if doing a GA calculation, return absolute value negative
            % (find most negative value), override other fit_values
            fit_value = -abs(fit_value);

            % if doing an eff calculation (i.e. actually doing GA), save the output
            if param.saveIntermediate
                temp = [vars output.MV output.MV2 output.MX output.MX2 output.flux_right output.flux_left output.en fit_value];
                % save intermediate data to a text file
                f = fullfile('defectOutput',['dynOpt_intermediate_' param.job_title '_' param.tor '.int']);
                save(f,'temp','-ascii','-append'); 
            end
        end
    else
        % return raw output
        fit_value = output;
    end

    %% toc
    totalEnd = toc(totalStart);
    fprintf('Finished a dynamics fitness run: %f seconds\n',totalEnd);
end