function [Psi,output] = propagate_1D_lvn(Psi,param)
%% Density matrix ODE solver

% note that the ode solver only accepts vector input.

% k-space transformation convention: fft columns, ifft the rows for the
% forward transform. ifft columns, fft rows for the inverse transform. The
% actual order is not important as pairs of fourier transforms are norm
% preserving

    % null output
    psi_pop    = [];
    MV         = [];
    MV2        = [];
    MX         = [];
    MX2        = [];
    flux_right = [];
    flux_left  = [];
    en         = [];
    % estimate on how long the data will be
    Nt_est = ceil((param.simul_time*(1-param.cutoff))/param.dt);
    if param.doNorm
        % <Psi|Psi>
        psi_pop = zeros(Nt_est,4);
        if param.doNorm == 2
            [~,imax] = max(sum(param.V,2));
            [~,imin] = min(sum(param.V,2));
            % determine easy and hard directions
            easy = param.X<param.X(imax) | param.X>param.X(imin);
            hard = param.X>param.X(imax) & param.X<param.X(imin);
        end
    end
    if param.doMV
        % <p>/m
        MV = zeros(Nt_est,1);
    end
    if param.doMV2
        % <p^2>/m^2
        MV2 = zeros(Nt_est,1);
    end
    if param.doX
        % <x>
        MX = zeros(Nt_est,1);
    end
    if param.doX2
        % <x^2>
        MX2 = zeros(Nt_est,1);
    end
    if param.doFlux
        % flux, counted as prob density on that matrix element
        flux_right = zeros(Nt_est,1);
        flux_left = zeros(Nt_est,1);
    end
    if param.doEnergy == 1
        % <T> only
        en = zeros(Nt_est,1);
    elseif param.doEnergy == 2
        % <T> and <V>
        en = zeros(Nt_est,2);
    end
    time_vector = zeros(Nt_est,1);
    i_ = 1; % counter for output
    
    % k-representation of velocity operators
    opV  = param.kxs/param.effMass; % MV, flux (derivative)
    opV2 = param.kxs.^2/param.effMass^2; % MV2 (2nd derivative)
    
    % coordinate operators
    opX  = param.X;
    opX2 = param.X.^2;
    
    % precalculate the kinetic difference matrix
    % [T,p]_ij = (T_i - T_j)p_ij
    Tcomm = bsxfun(@minus,param.T,param.T.');
    
    % preallocate the potential
    V = sum(param.V,2); % if ratchet is off, use fully on state;
    Vcomm = bsxfun(@minus,V,V.');
    
    if param.doBath(1) >= 2
        % preallocate lindblad operators
        
        % [x,rho]
        Xcomm = bsxfun(@minus,param.X,param.X.');
        % correct Xcomm by periodic boundary condition
        Xcomm(Xcomm>param.Lx/2) = Xcomm(Xcomm>param.Lx/2) - param.Lx;
        Xcomm(Xcomm<-param.Lx/2) = Xcomm(Xcomm<-param.Lx/2) + param.Lx;
        
        % [p,rho]
        Pcomm = bsxfun(@minus,param.kxs,param.kxs.');
        
        % prefactor [x,{p,rho}] iy/h
        % {p,rho}
        Panticomm = bsxfun(@plus,param.kxs,param.kxs.');
        xprho = 1i*param.gamma*Xcomm;
        
        % prefactor [x,[x,rho]] 2mykT/h^2
        xxrho = 2*param.effMass*param.gamma*param.kT*Xcomm.^2;
        
        % prefacotr [p,[p,rho]] y/8mkT
        pprho = param.gamma/(8*param.effMass*param.kT)*Pcomm.^2;
    end
    
    % RHS of ODE
    function dpdt = lvn_RHS(t,p)
        
        if param.doRatchet
            % grab the timing information
            Vt = propagate_1D_Vt(t,param);

            % form the full potential (column vectors of V * column vector vt)
            V = param.V*Vt;
            
            % calculate the potential difference matrix, vectorized
            Vcomm = bsxfun(@minus,V,V.');
        end
        
        % need to turn p back into a matrix for fft to work
        p = reshape(p,param.Nx,param.Nx);
        
        % move p into k-space
        kp = ifft(fft(p),[],2);
        
        % Actually calculate the liouvillian
        switch param.doBath(1)
            case 0
                % unitary evolution
                dpdt = -1i*(fft(ifft(Tcomm.*kp),[],2) + Vcomm.*p);
            case 1
                % (S)SH
            case 2
                % Lindblad - quadratic free particle
                dpdt = -1i*(fft(ifft(Tcomm.*kp),[],2) + Vcomm.*p) ... % unitary evolution
                       - xprho.*(fft(ifft(Panticomm.*kp),[],2)) ... % i\h y [x,{p,rho}]
                       - xxrho.*p... % 2mykT/h^2[x,[x,rho]]
                       - fft(ifft(pprho.*kp),[],2); % y/8mkT[p,[p,rho]]
            case 3
                % Lindblad - decoherence only term (not sure if correct)
                dpdt = -1i*(fft(ifft(Tcomm.*kp),[],2) + Vcomm.*p) ...
                       - xxrho.*p; % 2mykT/h^2[x,[x,rho]]
        end
        
        % turn it back into a vector for ode45
        dpdt = dpdt(:);
    end

    % output function called at EVERY successful integration
    % collects observables
    function status = myOutputFcn(t,p,flag)
        status = 0; % continue the integration (no stop conditions here)
        
        % <O> = tr(pO) = tr(Op)
        % If O is diagonal, it is simply element by element multiplication
        % of O and the diagonal of p. However, the velocity operator has to
        % actually be properly multiplied in k-space since the fourier
        % transform is non-local and the result is only normalized in
        % fft/ifft pairs

        % cut off time to discard earlier simulation time
        if t >= param.simul_time*param.cutoff + param.continue
            if strcmp(flag,'init')
                % t is given as tspan so take the first element only
                t = t(1);
            end
            if isempty(flag) || strcmp(flag,'init')
                % calculate observables at time = t

                % reshape from vector to matrix
                p = reshape(p,param.Nx,param.Nx);

                % move p into k-space
                kp = ifft(fft(p),[],2);

                % all ops in k-space are diag, so only need .*
                kp = diag(kp);

                % store the current time
                time_vector(i_) = t;

                if param.doMV
                    % <p>/m (diag turns back into square for i/fft)
                    MV(i_) = trace(fft(ifft(diag(opV.*kp)),[],2));
                end

                if param.doMV2
                    % <p^2>/m^2
                    MV2(i_) = trace(fft(ifft(diag(opV2.*kp)),[],2));
                end

                if param.doX
                    % <x>
                    MX(i_) = sum(diag(p).*opX);
                end

                if param.doX2
                    % <x^2>
                    MX2(i_) = sum(diag(p).*opX2);
                end

                if param.doFlux
                    % flux, counted as prob density on that matrix element
                    flux_right(i_) = p(fPoint+1,fPoint+1);
                    flux_left(i_)  = p(param.Nx-fPoint,param.Nx-fPoint);
                end

                if param.doEnergy
                    % calculate the energy at every time step  
                    % T
                    if param.doMV2
                        % <T> = <p^2/2m> = <p^2>/2m = <v^2>*m/2
                        en(i_,1) = MV2(i_)*param.effMass/2;
                    else
                        en(i_,1) = trace(fft(ifft(diag(param.T.*kp)),[],2));
                    end
                    if param.doEnergy == 2
                        en(i_,2) = sum(diag(p).*V);
                    end
                end
                
                if param.doNorm
                    % <Psi|Psi>
                    psi_pop(i_,1) = trace(p); % prob density
                    psi_pop(i_,2) = sum(p(:)) - psi_pop(i_,1); % off diag
                    if param.doNorm == 2
                        % find <Psi|Psi> with energy below V
                        ind = en(i_) < V;
                        % easy direction (less length to peak)
                        psi_pop(i_,3) = sum(p(diag(ind & easy)));
                        % hard direction (more distance to peak)
                        psi_pop(i_,4) = sum(p(diag(ind & hard)));
                    end
                end

                % increment observables counter
                i_ = i_ + 1;

                if t == param.simul_time + param.continue
                    % return the final density matrix
                    Psi = p;
                end
            end
        end
        if strcmp(flag, 'done')
            % discard unused space
            psi_pop(i_:end,:)   = [];
            MV(i_:end)          = [];
            MV2(i_:end)         = [];
            MX(i_:end)          = [];
            MX2(i_:end)         = [];
            flux_right(i_:end)  = [];
            flux_left(i_:end)   = [];
            en(i_:end,:)        = [];
            time_vector(i_:end) = [];
            
            % condition outputs to be real
            output.psi_pop     = psi_pop; % don't take real to see if propagation is acting funny
            output.MV          = real(MV);
            output.MV2         = real(MV2);
            output.MX          = real(MX);
            output.MX2         = real(MX2);
            output.flux_right  = flux_right;
            output.flux_left   = flux_left;
            output.en          = real(en);
            output.time_vector = time_vector;
            
            if param.doMovie(1)
                % clear all generated figure thingees
                close(f)
            end
        end
    end

    % event function for making movies
    function [value,isterminal,direction] = myEventsFcn(t,p)
        % integration does not stop in events function, it only captures
        % times when the movie should be made
        isterminal = 0;
        
        % uses a sine to locate zeros
        direction = 0;
        
        % use a sine with periodicity of modFrames
        value = sin(2*pi*t/tauFrame);
        
        if t >= param.simul_time*(1-param.doMovie(1)) + param.continue
            % just handle the movie making in the events function
            if abs(value) <= 5E-12; % can't get to zero
                % reshape from vector to matrix
                p = reshape(p,param.Nx,param.Nx);

                % make the movie
                if param.doMovie(2) == 1
                    h.CData = real(p);
                    caxis([0 max(real(p(:)))]);
                elseif param.doMovie(2) == 2
                    if param.doRatchet
                        set(h(1),'YData',V);
                    end
                    h(2).YData = real(diag(p));
                    th.String = t;
                end
                
                png_name = ['movie_imgs\timeinterval_' sprintf('%05d',m_)];
                print(f,'-dpng', png_name);
                m_ = m_ + 1;
            end
        end
    end

    %% build the density matrix |Psi><Psi|
    if size(Psi,1) == size(Psi,2)
        % if the Psi provided is a square matrix, assume it is a density
        % matrix and don't overwrite it with the fn's internal starting
        % conditions
        p0 = Psi;
    else
        if param.bias < 0; 
            % indicates a thermal state is needed
            switch fix(param.bias)
                case {-1,-4}
                    KE = param.T(1:param.Nx/2+1);
                case {-2,-5} % left+right traveling
                    KE = param.T;
                case {-3,-6}
                    KE = [0 param.T(param.Nx/2+2:end)];
                case -9
                    KE = Psi(:,end);
            end

            % partition, kT = elecKE 
            Z = exp(-KE/param.elecKE);

            p0 = zeros(param.Nx);
            for a = 1:numel(KE)
                if fix(param.bias) == -9
                    % thermal eigenfunction wavefunction
                    p0 = p0 + Z(a)*(Psi(:,a)*Psi(:,a)');
                else
                    % create the wavefunction at the particular kinetic energy
                    Psi = initializeWF_1D(param,KE(a));
                    % add it to the density matrix with the right probability
                    p0 = p0 + Z(a)*(Psi*Psi');
                end
            end
            p0 = p0/sum(Z); % normalize the weightings by Tr(exp(-bH))
        else
            % pure state density matrix
            p0 = Psi*Psi';
        end

        % condition the diagonal to be real
        p0(1:param.Nx+1:end) = real(diag(p0));
    end
    
    %% options
    % options for ode45
    options = odeset;
    options = odeset(options,'MaxStep',param.dt); % maximum step is highest freq of system
    options = odeset(options,'InitialStep',param.dt); % try initial step of dt
    options = odeset(options,'OutputFcn',@myOutputFcn);
    options = odeset(options,'Refine',1); % ode45 has default at 4
    options = odeset(options,'Stats','on');
    if param.doMovie(1)
        % enable the event location function
        fprintf('movie!')
        m_ = 1;
        
        % enable events function
        options = odeset(options,'Events',@myEventsFcn);
        
        % fraction of Nt to visualize, 1-param.doMovie
        frames = 1000;
        tauFrame = round(param.simul_time*param.doMovie(1)/frames);
        
        f = figure(1);
        clf(f);
        f.Visible = 'off';
        f.PaperUnits = 'inches';
        f.PaperPosition = [0 0 4 3];
        
        if param.doMovie(2) == 1
            fprintf(' - full density matrix\n')
            % set up the density matrix imagesc
            h = imagesc(p0);
            set(gca,'YDir','normal');
            axis square;
            figwu;
        elseif param.doMovie(2) == 2  
            fprintf(' - diag density matrix\n')
            % Set up a 2 row matrix
            plot_data(1:param.Nx,1) = sum(param.V,2);
            plot_data(1:param.Nx,2) = diag(p0);
            % make a 2 plot object
            maxV = max(param.V(:));
            minV = min(param.V(:));
            h = plot(param.X*param.a0/1E-9,plot_data); % in nm
            ylim([minV max([maxV max(diag(p0)) 0.08])]); 
                % 0.08 magic number for highest of rho during simulation
            xlabel('Length (a0)');
            % set the colour
            h(1).Color = 'm';
            h(2).Color = 'b';
            % make a handle to the title
            th = title('t = 0','FontWeight','normal');
            % zero line
            line([param.X(1) param.X(end)]*param.a0/1E-9,[0 0],'LineStyle','--');
        end
    end
    
    %% propagate the density matrix with ode45
    % don't need results or time as the output function will grab values during propagation
    if size(Psi,1) == size(Psi,2)
        % continue from where the simulation left off
        prop_time = [param.continue param.continue+param.simul_time];
    else
        % fresh propagation
        prop_time = [0 param.simul_time];
    end
    switch param.ode_solver
        case 23
            ode23(@lvn_RHS,prop_time,p0,options);
        case 45
            ode45(@lvn_RHS,prop_time,p0,options);
        case 113
            ode113(@lvn_RHS,prop_time,p0,options);
    end
end