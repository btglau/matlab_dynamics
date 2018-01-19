function [ param ] = initializeTime_1D(param)
%% Initialize time parts of the simulation

    % set the time step to the limit of the grid
    param.dt = floor(100/(max(abs(param.V(:))) + pi^2/(2*param.effMass*param.dx^2)))/100;
    
    % determine the total time
    param.all_taus = zeros(param.unique_tiles,1);
    for a = 1:param.unique_tiles
        ms = param.time_matrix(param.time_ind(a):param.time_ind(a+1)-1);
        switch ms(1)
            case -1
                param.all_taus(a) = sum(ms(2:end-1));
            case -2
                param.all_taus(a) = ms(2);
        end
    end
    
    if param.doRatchet 
        % if doing the ratchet, need to set some minimum time to allow e-
        % to reach a steady state
        if param.doBath(1) == 0
            % 5 ps or x times the max period, whichever is greater
            param.simul_time = max([2.0678E5 param.time_multiplier*max(param.all_taus)]);
        elseif param.doBath(1) > 0
            % if doing a bath, max 2 ps or x times period
            param.simul_time = max([1.0339E5 param.time_multiplier*max(param.all_taus)]);
            if param.simul_time > 2.0678E5
                % if simul time is greater than 5 ps, drop the multiplier
                param.simul_time = 10*max(param.all_taus);
            end
        end
    else
        % if doing a static potential, ignore minimum time
        param.simul_time = param.time_multiplier*max(param.all_taus);
    end
    
    % set the number of time steps
    param.Nt = ceil(param.simul_time/param.dt);
    
    % get the maximum of any time fourier series so it can be normalized
    if any(param.time_matrix == -2)
        % build the fourier series (yeah .. not most elegant ..)
        param.fourier_norm = ones(sum(param.time_matrix == -2),1);
        i_ = 1;
        for a = 1:param.unique_tiles
            ms = param.time_matrix(param.time_ind(a):param.time_ind(a+1)-1);
            switch ms(1)
                case -1
                    % don't do anything
                case -2
                    % chop out the -2 and tau to make indexing in loop easier
                    ms = ms(3:end);
                    
                    % if there are more than 1 fourier terms
                    if length(ms) > 2
                        % check over one period of 0:pi (since sin.^2)
                        per_t = linspace(0,pi,2^11).';

                        % fourier series of sin^2 terms
                        vt = sum(bsxfun(@times,ms(1:2:end),sin(bsxfun(@plus,bsxfun(@times,(1:length(ms)/2),per_t),ms(2:2:end))).^2),2);

                        % grab the maximum
                        param.fourier_norm(i_) = max(vt);

                        % counter that increments only for fourier time series
                        i_ = i_ + 1;
                    end
            end
        end
    end
end