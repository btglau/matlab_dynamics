function [Psi,output] = propagate_1D(Psi,param)
%% propagate a grid-based wavefunction with SPO, SIL + SH bath, or CHSH + SH bath

    % null output
    psi_pop     = zeros(0,1);
    MV          = zeros(0,1);
    MV2         = zeros(0,1);
    MX          = zeros(0,1);
    MX2         = zeros(0,1);
    flux_right  = zeros(0,1);
    flux_left   = zeros(0,1);
    en          = zeros(0,1);
    
    % precalculated values
    halfdt      = -1i*param.dt/2; % potential
    V           = sum(param.V,2); % if ratchet is off, use fully on state
    
    % Nt_cutoff: when to start collecting observables
    Nt_cutoff = round(param.Nt*param.cutoff);
    i_ = 1;
    
    if param.doRatchet
        % grab the timing information
        Vt = propagate_1D_Vt((0:param.Nt).'*param.dt,param);
    end
    
    switch param.prop
        % set up propagator parameters
        case 1
            expT = exp(-1i*param.T*param.dt); % full time step kinetic operator
            if param.doNIP
                expV = exp((V  + param.imagV)*halfdt);
            else
                expV = exp(V * halfdt);
            end
        case 2
            % Perform short iteractive lanczos propagation method
            order = 5;
            krylov_h_subspace = zeros(order);
            % squeeze to overload for bath on / off
            krylov_basis = squeeze(zeros([size(Psi) order]));
        case 3
            % no need to do anything
    end
    
    if param.doNorm
        % <Psi|Psi>
        psi_pop = zeros(param.Nt-Nt_cutoff,1);
    end
    
    if param.doMV || param.doMV2
        if param.doMV
            % <p>/m
            opV = param.kxs/param.effMass; % MV, flux (derivative)
            MV = zeros(param.Nt-Nt_cutoff,1);
        end
        if param.doMV2
            % <p^2>/m^2
            opV2 = param.kxs.^2/param.effMass^2; % MV2 (2nd derivative)
            MV2 = zeros(param.Nt-Nt_cutoff,1);
        end
    end

    if param.doX
        % <x>
        opX = param.X;
        MX = zeros(param.Nt-Nt_cutoff,1);
    end
    
    if param.doX2
        % <x^2>
        opX2 = param.X.^2;
        MX2 = zeros(param.Nt-Nt_cutoff,1);
    end
    
    if param.doFlux
        % do flux UNFUDGED (need to divide sum of flux by dx)
        flux_right = zeros(param.Nt-Nt_cutoff,1);
        flux_left = zeros(param.Nt-Nt_cutoff,1);
    end
    
    if param.doEnergy > 0
        % calculate the energy at every time step
        en = zeros(param.Nt-Nt_cutoff,2);
    end
    
    if param.doMovie
        
        fprintf('movie!\n')
        
        maxV = max(param.V(:));
        minV = min(param.V(:));
        
        f = figure('Visible','off');
        set(f,'PaperUnits','inches','PaperPosition',[0 0 6 4]);
        % Set up a 2 row matrix
        plot_data(1:param.Nx,1) = V;
        plot_data(1:param.Nx,2) = sum(abs(Psi).^2,2);
        % make a 2 plot object
        h = plot(param.X,plot_data);
        ylim([minV maxV]);
        xlabel('Length (a0)');
        % set the colour
        set(h(1),'Color','m');
        set(h(2),'Color','b');
        y_ = 1;
        
        % fraction of Nt to visualize
        numFrame = round(param.Nt*param.doMovie);
        modFrame = round(param.Nt*param.doMovie/1000);
    end
    
    % Check for GPUs!!!
    if gpuDeviceCount > 0
        %g = gpuDevice(); % get the default GPU device
        %Psi = gpuArray(Psi);
        %fprintf(['Connected to: ' g.Name ', for GPU calculation!\n'])
    else
        %fprintf('No GPU detected!')
    end

    % N propagations ------------------------------------------------------
    for a = 1:param.Nt
        
        if a > Nt_cutoff % only collect data past the cutoff
            % Precalculate for multiple observable calculations
            conjPsi = conj(Psi);

            if param.doNIP
                % if absorbing potential, renormalize the wavefunction
                % (commented out because author feels it is physical to look
                % only at the remaining probability density)
                %norm_const = 1/sqrt(sum(abs(Psi(:)).^2));
                %conjPsi = conjPsi*norm_const;
            end

            if param.doNorm
                % <Psi|Psi>
                psi_pop(i_) = sum(conjPsi(:).*Psi(:));
            end

            if param.doMV || param.doMV2
                % put Psi into k-space
                kPsi = fft(Psi);

                if param.doMV
                    % <p>/m
                    if param.doBath(1) == 0
                        PsiOpPsi = conjPsi.*ifft(opV.*kPsi);
                    else
                        PsiOpPsi = conjPsi.*ifft(bsxfun(@times,opV,kPsi));
                    end
                    MV(i_) = sum(PsiOpPsi(:));
                end
                if param.doMV2
                    % <p^2>/m^2
                    if param.doBath(1) == 0
                        PsiOpPsi = conjPsi.*ifft(opV2.*kPsi);
                    else
                        PsiOpPsi = conjPsi.*ifft(bsxfun(@times,opV2,kPsi));
                    end
                    MV2(i_) = sum(PsiOpPsi(:));
                end
            end

            if param.doX
                % <x>
                if param.doBath(1) == 0
                    PsiOpPsi = conjPsi.*opX.*Psi;
                else
                    PsiOpPsi = conjPsi.*bsxfun(@times,opX,Psi);
                end
                MX(i_) = sum(PsiOpPsi(:));
            end

            if param.doX2
                % <x^2>
                if param.doBath(1) == 0
                    PsiOpPsi = conjPsi.*opX2.*Psi;
                else
                    PsiOpPsi = conjPsi.*bsxfun(@times,opX2,Psi);
                end
                MX2(i_) = sum(PsiOpPsi(:));
            end

            if param.doFlux
                % do flux UNFUDGED
                if param.doBath(1) == 0
                    kPsi        = ifft(1i*param.kxs.*fft(Psi));
                else
                    kPsi        = ifft(bsxfun(@times,1i*param.kxs,fft(Psi)));
                end
                flux_right(i_)   = sum(imag(conj(Psi(param.fPoint+1,:)).*kPsi(param.fPoint+1,:)));
                flux_left(i_)    = sum(imag(conj(Psi(param.Nx-param.fPoint,:)).*kPsi(param.Nx-param.fPoint,:)));
            end

            if param.doEnergy
                % calculate the energy at every time step
                % T
                if param.doBath(1)
                    PsiOpPsi = conjPsi.*ifft(bsxfun(@times,param.T,fft(Psi)));
                else
                    PsiOpPsi = conjPsi.*ifft(param.T.*fft(Psi));
                end
                en(i_,1) = sum(PsiOpPsi(:));
                if param.doEnergy ==2
                    % V
                    if param.doBath(1)
                        PsiOpPsi = conjPsi.*bsxfun(@times,V,Psi);
                    else
                        PsiOpPsi = conjPsi.*V.*Psi;
                    end
                    en(i_,2) = sum(PsiOpPsi(:));
                end
            end
            i_ = i_ + 1;
        end
        
        if param.doMovie
            % a dummy counter for the movie frame count
            % update h(2) to the current time step and print the png
            if mod(a,modFrame) == 0 && a <= numFrame
                if param.doRatchet
                    set(h(1),'YData',V);
                end
                set(h(2),'YData',sum(abs(Psi).^2,2));
                png_name = ['movie_imgs\timeinterval_' sprintf('%05d',y_)];
                print('-dpng', png_name);
                y_ = y_ + 1;
            end
        end

        if param.doRatchet
            % update the time dependent potential
            V = param.V*Vt(a);
        end
        
        switch param.prop
            case 1
                if param.doBath
                    error('SPO cannot propagate with a bath')
                end
                    % SPO only if the bath is off
                    % form the new potential
                    if param.doRatchet
                        if param.doNIP
                            expV = exp((V  + param.imagV)*halfdt);
                        else
                            expV = exp(V * halfdt);
                        end
                    end

                    Psi = expV.*Psi; % potential
                    Psi = ifft(expT.*fft(Psi));
                    Psi = expV.*Psi;
            case 2
                if param.doBath(1) == 0
                    % first basis function q_0
                    krylov_basis(:,1) = Psi;

                    % second basis function q_1 (iteration b = 2)
                    Hqi = ifft(param.T.*fft(krylov_basis(:,1))) + V.*krylov_basis(:,1); % Hq_0
                    krylov_h_subspace(1,1) = sum(conj(krylov_basis(:,1)).*Hqi); % alpha_0
                    krylov_basis(:,2) = Hqi - krylov_h_subspace(1,1)*krylov_basis(:,1); % Hq_0 - alpha_0 q_0
                    if 1 == 1
                        norm_const = 1/norm(krylov_basis(:,2));
                        krylov_h_subspace(1,2) = sum(conj(krylov_basis(:,2)*norm_const).*Hqi); % use normalized q_1 to find beta_0
                        krylov_basis(:,2) = krylov_basis(:,2)*(1/krylov_h_subspace(1,2)); % q1 = () / beta_0
                    else
                        krylov_basis(:,2) = krylov_basis(:,2)/norm(krylov_basis(:,2),2); % normalize q_1
                        krylov_h_subspace(1,2) = sum(conj(krylov_basis(:,2).*Hqi)); % beta_0
                    end

                    % loop through the rest of the basis up to order (=b-1)
                    for b = 3:order
                        Hqi = ifft(param.T.*fft(krylov_basis(:,b-1))) + V.*krylov_basis(:,b-1); % Hq_(b-2)
                        krylov_h_subspace(b-1,b-1) = sum(conj(krylov_basis(:,b-1)).*Hqi); % alpha_(b-2)
                        krylov_basis(:,b) = Hqi - krylov_h_subspace(b-2,b-1)*krylov_basis(:,b-2) - krylov_h_subspace(b-1,b-1)*krylov_basis(:,b-1); % Hq_(b-2) - beta_(b-3)q_(b-3) - alpha(b-2)q_(b-2)
                        if 1 == 1
                            norm_const = 1/norm(krylov_basis(:,b));
                            krylov_h_subspace(b-1,b) = sum(conj(krylov_basis(:,b)*norm_const).*Hqi); % beta_{b-2} with normalized q_{b-1}
                            krylov_basis(:,b) = krylov_basis(:,b)*(1/krylov_h_subspace(b-1,b)); % q_{b-1} = () / beta_{b-2}
                        else
                            krylov_basis(:,b) = krylov_basis(:,b)/norm(krylov_basis(:,b),2); % normalize q_{b-1}
                            krylov_h_subspace(b-1,b) = sum(conj(krylov_basis(:,b).*Hqi)); % beta_{b-2}
                        end
                    end

                    % fill in the last part of the krylov hamiltonian
                    Hqi = ifft(param.T.*fft(krylov_basis(:,b))) + V.*krylov_basis(:,b);
                    krylov_h_subspace(order,order) = sum(conj(krylov_basis(:,b)).*Hqi);

                    % reflect the krylov hamiltonian (hermitian)
                    krylov_h_subspace(2:order+1:end) = krylov_h_subspace(order+1:order+1:end);
                    
                    % diagonalize the krylov hamiltonian
                    % Z'*krylov_h_subspace*Z = D
                    [Z,D] = eig(krylov_h_subspace);

                    % construct the propagator in the krylov subspace
                    % only the first column contributes to time propagation
                    Ut = Z * diag(exp(-1i*diag(D)*param.dt)) * Z';

                    % calculate the new propagated wavefunction Psi
                    Psi = sum(bsxfun(@times,Ut(:,1).',krylov_basis),2);
                elseif param.doBath(1) == 1
                    % first basis function q_0
                    krylov_basis(:,:,1) = Psi;

                    % second basis function q_1 (iteration b = 2)
                    Hqi = propagate_1D_HSB(param,param.T,V,krylov_basis(:,:,1)); % Hq_1
                    krylov_h_subspace(1,1) = sum(sum(conj(krylov_basis(:,:,1)).*Hqi)); % alpha_0
                    krylov_basis(:,:,2) = Hqi - krylov_h_subspace(1,1)*krylov_basis(:,:,1); % unnormalized krylov basis vector q_1
                    norm_const = 1/sqrt(sum(sum(abs(krylov_basis(:,:,2)).^2)));
                    krylov_h_subspace(1,2) = sum(conj(krylov_basis(:,:,2)*norm_const).*Hqi); % use normalized q_1 to find beta_0
                    krylov_basis(:,:,2) = krylov_basis(:,:,2)*(1/krylov_h_subspace(1,2)); % q1 = () / beta_0

                    % loop through the rest of the basis up to order (=b-1)
                    for b = 3:order
                        Hqi = propagate_1D_HSB(param,param.T,V,krylov_basis(:,:,b-1)); % Hq_(b-2)
                        krylov_h_subspace(b-1,b-1) = sum(sum(conj(krylov_basis(:,:,b-1)).*Hqi)); % alpha_{b-2}
                        krylov_basis(:,:,b) = Hqi - krylov_h_subspace(b-2,b-1)*krylov_basis(:,:,b-2) - krylov_h_subspace(b-1,b-1)*krylov_basis(:,:,b-1); % unnormalized krylov basis vector {b-1}
                        norm_const = 1/sqrt(sum(sum(abs(krylov_basis(:,:,b)).^2)));
                        krylov_h_subspace(b-1,b) = sum(conj(krylov_basis(:,:,2)*norm_const).*Hqi); % beta_{b-2}
                        krylov_basis(:,:,b) = krylov_basis(:,:,b)*(1/krylov_h_subspace(b-1,b)); % q_{b-1} = () / beta_{b-2}
                    end

                    % fill in the last part of the krylov hamiltonian
                    Hqi = propagate_1D_HSB(param,param.T,V,krylov_basis(:,:,b));
                    krylov_h_subspace(order,order) = sum(sum(conj(krylov_basis(:,:,b)).*Hqi));

                    % reflect the krylov hamiltonian
                    krylov_h_subspace(2:order+1:end) = krylov_h_subspace(order+1:order+1:end);

                    % diagonalize the krylov hamiltonian
                    % Z'*krylov_h_subspace*Z = D
                    [Z,D] = eig(krylov_h_subspace);

                    % construct the propagator in the krylov subspace
                    % only the first column contributes to time propagation
                    Ut = Z * diag(exp(-1i*diag(D)*param.dt)) * Z';

                    % calculate the new propagated wavefunction Psi
                    Psi = sum(bsxfun(@times,reshape(Ut(:,1),[1 1 order]),krylov_basis),3);
                end
            case 3
                if param.doRatchet
                    % no need to update variables if ratchet is off
                    Vmin    = min(V);
                    dEgrid  = (pi/param.dx)^2/(2*param.effMass) + max(V) - Vmin;
                    chshT   = 2*param.T/dEgrid;
                    chshV   = V - dEgrid/2 - Vmin;
                    chshV   = 2*chshV/dEgrid;
                    alpha   = dEgrid*param.dt/2;
                    
                    % grab the bessel functions
                    [chshBessel,N] = propagate_chsh_bessel(alpha);
                end
                
                if param.doBath(1) == 0
                    % phi_0 is the original wavefunction 'psi'
                    phi_0 = Psi;

                    % apply the hamiltonian to phi_0 to generate phi_1
                    phi_1 = -1i*(ifft(chshT.*fft(phi_0)) + chshV.*phi_0);

                    % phi_2 is generated in part by applying the hamiltonian to phi_1
                    phi_2 = -2i*(ifft(chshT.*fft(phi_1)) + chshV.*phi_1) + phi_0;
                    
                    % Start summing the chsh expansion
                    Psi = chshBessel(3)*phi_2 + chshBessel(2)*phi_1 + chshBessel(1)*phi_0;

                    % Now loop the recursion
                    for d = 4:N % this index is for the bessel functions

                        phi_0 = phi_1; % phi_{n-1}
                        phi_1 = phi_2; % phi_{n}

                        % apply the hamiltonian to get phi_{n+1}
                        phi_2 = -2i*(ifft(chshT.*fft(phi_1)) + chshV.*phi_1) + phi_0;
                        Psi = Psi + chshBessel(d)*phi_2;
                    end
                elseif param.doBath(1) == 1
                    % phi_0 is the original wavefunction 'psi'
                    phi_0 = Psi;

                    % apply the hamiltonian to phi_0 to generate phi_1
                    phi_1 = -1i*propagate_1D_HSB(param,chshT,chshV,phi_0);

                    % phi_2 is generated in part by applying the hamiltonian to phi_1
                    phi_2 = -2i*propagate_1D_HSB(param,chshT,chshV,phi_1) + phi_0;
                    
                    % Start summing the chsh expansion
                    Psi = chshBessel(3)*phi_2 + chshBessel(2)*phi_1 + chshBessel(1)*phi_0;

                    % Now loop the recursion
                    for d = 4:N % this index is for the bessel functions

                        phi_0 = phi_1; % phi_{n-1}
                        phi_1 = phi_2; % phi_{n}      

                        % apply the hamiltonian to get phi_{n+1}
                        phi_2 = -2i*propagate_1D_HSB(param,chshT,chshV,phi_1) + phi_0;
                        Psi = Psi + chshBessel(d)*phi_2;
                    end
                end
                
                % phase correction
                Psi = Psi*exp(-1i*(Vmin*param.dt + alpha));
        end
        % end of a single time step
    end
    
    if gpuDeviceCount > 0
        %reset(g); % reset the gpuDevice
        %Psi = gather(Psi);
        %fprintf(['Finished calcs and resetting GPU: ' g.Name '\n']);
    end
    
    % close all figure frames to free up figure memory
    if param.doMovie;close all;end
    
    % create the output file (flux is fudged here)
    output.psi_pop      = psi_pop;
    output.MV           = real(MV);
    output.MV2          = real(MV2);
    output.MX           = real(MX);
    output.MX2          = real(MX2);
    output.flux_right   = flux_right/param.dx;
    output.flux_left    = flux_left/param.dx;
    output.en           = real(en);
end