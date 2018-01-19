function fit_value = mcFitness(vars,param)
% this function is a propagator + observables calculator
% which uses an analytic Newtonian expression to calculate electron
% trajectories, and then uses a monte carlo algorithm to select scattering 
% mechanisms and states.
%
% design of the propagation: the user supplies a total simulation time, and
% how many (equally spaced) points along this time at which to gather
% observables. since electrons undergo random free flight, the electron's time
% advances stochastically according to the self-scattering algorithm
% therefore, electrons are propagated forever until they hit the simulation
% time specified, then deleted from the array. the while loop continues while
% the array is non zero in number of electrons.
%
% for this program to work, the data collection time intervals should be
% appreciably greater than the monte carlo time steps, or else there is a
% chance that the electron's time will tick over 2 or more data collection
% intervals.

fprintf('\n==========================\n')
fprintf('the "Manchego Carlo"!\n')
fprintf('==========================\n')
    
    % F1 and G1 integrands (for acoustic phonon scattering)
    function y = F1(x)
        y = x.^2 ./(exp(x)-1);
    end
    function y = G1(x)
        y = x.^2 .* (1./(exp(x)-1) + 1);
    end

    function G = pgrad(grad,q)
        % function that returns 0->1 (for impurity density, or composition gradient)
        switch grad(1)
            case -1
                % rectangles
                % number of little rectangles
                nr = numel(grad(3:end));
                % how wide each little rectangle is
                ru = grad(2)/nr;
                % since the gradient is periodic, mod everything by the 
                % period of one stepwise gradient
                q = mod(q,grad(2));
                % bin mod(q) into little rectangles, then assign values
                % from the gradient provided
                G = discretize(q,0:ru:ru*nr,grad(3:end));
            case -2
                % fourier series
                % imp(1) = 2pi / nL, n = repeat units
                grad = grad'; % transpose grad as q is a column vector
                insine = (1:numel(grad(3:end))/2) * grad(2); % fourier index = 1,2,3,...
                G = bsxfun(@times,insine,q); % 2pi/L * insine * q
                G = sin(bsxfun(@plus,grad(4:2:end),G)); % add the phase, then sine
                G = sum(bsxfun(@times,grad(3:2:end),G),2); % multiply coeffs, then sum
            case -3
                % simple polynomial evaluation
                q = mod(q,grad(2))/grad(2);
                G = polyval(grad(3:end),q);
            case -4
                % sawtooth
                q = mod(q,grad(2))/grad(2); % from 0 -> 1
                q_ind = q < grad(3);
                q(q_ind) = q(q_ind)/grad(3); % x/a
                q(~q_ind) = (q(~q_ind)-1)/(grad(3)-1); %  (x-1)/(a-1)
                G = q;
            otherwise
                % imp(1) = 0, no scattering profile
                q(:) = 0;
                G = q;
        end
    end

    function [grad,g_norm,g_min,g] = grad_normalize(grad,L)
        % normalize a gradient (composition or impurity)
        % grad: [x;y;z]
        % L: [Lx Ly Lz]
        g_min = 0;
        g_norm = 1;
        
        % convert the integer number of periods to distance per period
        % (useful for mod operations)
        for a_ = 1:size(grad,1)
            switch grad(1,a_)
                case {-1,-4}
                    % blocks or sawtooth
                    grad(2,a_) = L(a_)/grad(2,a_);
                case -2
                    % fourier series
                    grad(2,a_) = 2*pi/(L(a_)/grad(2,a_));
                case -3
                    % polynomial
                    grad(2,a_) = L(a_)/grad(2,a_);
                    g_endpoints = [polyval(grad(a_,3:end),grad(2,a_)) polyval(grad(a_,3:end),0)];
                    g_min = min(g_endpoints);
                    g_norm = 1/(max(g_endpoints)-g_min);
                % sawtooth (case -4) is automatically normalized to 1
            end
        end
        
        % more than one gradient, or fourier series
        if nnz(grad(1,:)) > 1 || any(grad(1,:)) == -2
            % find the max / min of the gradient through iterative grids
            grid_points = 128;
            grid_max = 512/nnz(grad(:,1));
            refine = 16; % additive refinement (not multiplicative)
            g_max = 1;
            g_max_prev = 0;
            g_min = 1;
            g_min_prev = 0;
            while abs((g_max-g_max_prev)/g_max_prev) > 1E-3 || abs((g_min-g_min_prev)/g_min_prev) > 1E-3
                % form the grid and the grid representation of the scattering profile
                [gx,gy,gz] = meshgrid(pgrad(grad(:,1),linspace(0,L(1),grid_points).'),...
                                      pgrad(grad(:,2),linspace(0,L(2),grid_points).'),...
                                      pgrad(grad(:,3),linspace(0,L(3),grid_points).'));
                g_full = gx + gy + gz;
                % shuffle some values around
                g_max_prev = g_max;
                g_min_prev = g_min;
                g_max = max(g_full(:));
                g_min = min(g_full(:));
                % refine the grid
                grid_points = grid_points + refine;
                if grid_points > grid_max
                    fprintf('\nDid not converge up to %g^3 points!\n',grid_max)
                    break
                end
            end

            % range of values
            g_norm = 1/(g_max - g_min);
        end
        
        % return a profile in x, y, and z
        g(:,1) = pgrad(grad(:,1),linspace(0,L(1),2^10)');
        g(:,2) = pgrad(grad(:,2),linspace(0,L(2),2^10)');
        g(:,3) = pgrad(grad(:,3),linspace(0,L(3),2^10)');
    end

    function dq = dqdt(k0,t0,ff,DC,AC)
        % this function returns dqdt
        % constant field
        % r'(t) = k0/m + E/m (t-t0)
        dq = k0/m_star + DC.*ff/m_star;
        if AC(2) ~= 0
            % time dependent field
            % r'(t) = k0/m + E/mw[cos(wto) - cos(wt)]
            dq = dq + AC(1)/(m_star*AC(2))*(cos(AC(2)*t0) - cos(AC(2)*(t0+ff)));
        end
    end

    function c = faster_fzero(q0,k0,t0,DC,AC,L,a,b)        
        % this (vectorized) function uses newton's method to 
        % quickly find the zero for the well-behaved function fun - calls 
        % ff_q. a/b are the two endpoints, L is an offset to make sure the
        % function is always crossing 0
        a_old = a;
        b_old = b;
        
        for i_ = 1:5 % use bisection to collapse endpoints
            c = (a+b)/2;
            fun_c = ff_q(q0,k0,t0,c,DC,AC)-L;
            % collapse endpoints a / b based on sign
            fun_a = ff_q(q0,k0,t0,a,DC,AC)-L;
            ffz_ind = fun_a.*fun_c < 0;
            b(ffz_ind) = c(ffz_ind);
            a(~ffz_ind) = c(~ffz_ind);
        end
        
        % newton-rhapson method
        c_old = (a+b)/2;
        c = c_old - (ff_q(q0,k0,t0,c_old,DC,AC)-L)./dqdt(k0,t0,c_old,DC,AC);
        c_ind = abs(c-c_old) > 1E-4 * abs(c);
        while any(c_ind)
            c_old(c_ind) = c(c_ind);
            c(c_ind) = c_old(c_ind) - (ff_q(q0(c_ind),k0(c_ind),t0(c_ind),c_old(c_ind),DC(c_ind),AC)-L(c_ind))...
                                ./dqdt(k0(c_ind),t0(c_ind),c_old(c_ind),DC(c_ind),AC);
            c_ind = abs(c-c_old) > 1E-5 * abs(c);
        end

        if any(c < a_old) || any(c > b_old)
            disp(a_old)
            disp(c)
            disp(b_old)
            error('newton''s method getting way off, outside [a b]!')
        end
    end

    % scattering probability function
    function P = scatProb(C,I,e_ke,valley)
        P = zeros(numel(e_ke),45);
        
        % impurity - S is normalized impurity gradient from 0 to ni
        if unionized %1-3
            % cold enough that impurities don't ionize
            w = e_ke/EB;
            % elastic scattering
            P(:,1) = 35.2./sqrt(w) .* (1+exp(-50*w)).*(1+80.6*w+23.7*w.^2)./(1+41.3*w+133*w.^2) .* ((log(1+w)./w) - (1+w/2-w.^2/6)./(1+w).^3);
            P(:,1) = P(:,1) .* I * a_0 / m_star;
            % inelastic 1 -> 2 excitation
            s_ind = e_ke > Enn;
            if any(s_ind)
                Unn = e_ke(s_ind)/Enn; % E/Enn
                P(s_ind,2) = 2*n^2/x * 1./Unn .* (1 - exp(-rnn*Unn)).*(Ann*(log(Unn) + 1./(2*Unn)) + (Bnn - Ann*log(2*n^2/x)).*(1 - 1./Unn))*pi*a_0^2;
                if isscalar(I)
                    % single value of I - no gradients
                    P(s_ind,2) = P(s_ind,2) .* sqrt(2*e_ke(s_ind)/m_star) * I;
                else
                    P(s_ind,2) = P(s_ind,2) .* sqrt(2*e_ke(s_ind)/m_star) .* I(s_ind);
                end
            end
            % inelastic ionization
            s_ind = w * n^2 > 1;
            if any(s_ind)
                Un = w(s_ind) * n^2;    
                P(s_ind,3) = 2*n^2./Un .* (1-exp(-rn*Un)) .* (An*log(Un) + (Bn - An*log(2*n^2))*(1-1./Un).^2)*pi*a_0^2;
                if isscalar(I)
                    P(s_ind,3) = P(s_ind,3) .* sqrt(2*e_ke(s_ind)/m_star) * I;
                else
                    P(s_ind,3) = P(s_ind,3) .* sqrt(2*e_ke(s_ind)/m_star) .* I(s_ind);
                end
            end
        else %1, 2-3 are empty
            % decided between BH and CW by comparing b_max to the screening
            % length, P(:,2) and P(:,3) are zero
            if BH
                P(:,1) = P_imp_BH*I.*sqrt(e_ke)./(1 + 4*e_ke/eps_b);
            else
                P(:,1) = P_imp_CW * sqrt(e_ke);
            end
        end
        
        % acoustic phonon absorption - done by interpolation
        P(:,4) = P_ac_abs(e_ke);
        % emission
        s_ind = e_ke > e_u;
        if any(s_ind)
            P(s_ind,5) = P_ac_em(e_ke(s_ind));
        end
        
        % six intervalley phonon processes possible
        % g-type processes have no difference between initial and final
        % valleys (even if the conduction band is split by strain)
        for p_ = 1:3
            % phonon absorption
            P(:,p_+5) = P_op(p_)*Nop(p_)*sqrt(e_ke + wop(p_));
            % emission
            ph_en = e_ke - wop(p_);
            s_ind = ph_en > 0;
            if any(s_ind)
                P(s_ind,p_+8) = P_op(p_)*(Nop(p_)+1)*sqrt(ph_en(s_ind));
            end
            if all(photon)
                % photon absorption with phonon absorption
                P(:,p_+27) = W_op(p_)*Nop(p_)*sqrt(e_ke + wop(p_) + photon(1)).*(2*e_ke + wop(p_) + photon(1));
                if photon(1) > wop(p_) % photon absorption w/ phonon emission
                    P(:,p_+30) = W_op(p_)*(Nop(p_)+1)*sqrt(e_ke + photon(1) - wop(p_)).*(2*e_ke + photon(1) - wop(p_));
                end
            end
        end
        % f-type processes of form 4->2 or 2->4 (0.67x (% Ge)(eV))
        for p_ = 4:6
            % absorption
            ph_en = e_ke + wop(p_) - C;
            s_ind = ph_en > 0;
            if any(s_ind)
                P(s_ind,p_+8) = P_op(p_)*Nop(p_)*sqrt(ph_en(s_ind));
            end
            % emission
            ph_en = e_ke - wop(p_) - C;
            s_ind = ph_en > 0;
            if any(s_ind)
                P(s_ind,p_+11) = P_op(p_)*(Nop(p_)+1)*sqrt(ph_en(s_ind));
            end
            if all(photon)
                ph_en = photon(1) + wop(p_) - C;
                s_ind = ph_en > 0;
                if any(s_ind)
                    P(s_ind,p_+30) = W_op(p_)*Nop(p_)*sqrt(e_ke(s_ind) + ph_en(s_ind)).*(2*e_ke(s_ind) + ph_en(s_ind));
                end
                ph_en = photon(1) - wop(p_) - C;
                s_ind = ph_en > 0;
                if any(s_ind)
                    P(s_ind,p_+33) = W_op(p_)*(Nop(p_)+1)*sqrt(e_ke(s_ind) + ph_en(s_ind)).*(2*e_ke(s_ind) + ph_en(s_ind));
                end
            end
        end
        % 4->2 f scattering only goes into 2 final valleys, not 4
        P(valley,[12:17 34:39]) = P(valley,[12:17 34:39])/2;
        % 4->4 f-type with 2 final valleys and no energy diff
        % divide by 2 as only 2 final valleys 4->4
        if any(valley)
            for p_ = 4:6
                % absorption
                P(valley,p_+14) = P_op(p_)*Nop(p_)*sqrt(e_ke(valley) + wop(p_))/2;
                % emission
                ph_en = e_ke - wop(p_);
                s_ind = ph_en > 0 & valley; % further filter for only electrons with energy to emit
                if any(s_ind)
                    P(s_ind,p_+17) = P_op(p_)*(Nop(p_)+1)*sqrt(ph_en(s_ind))/2;
                end
                if all(photon)
                    P(valley,p_+36) = W_op(p_)*Nop(p_)*sqrt(e_ke(valley) + wop(p_) + photon(1)).*(2*e_ke(valley) + wop(p_) + photon(1))/2;
                    if photon(1) > wop(p_) % photon absorption w/ phonon emission
                        P(valley,p_+39) = W_op(p_)*(Nop(p_)+1)*sqrt(e_ke(valley) + photon(1) - wop(p_)).*(2*e_ke(valley) + photon(1) - wop(p_))/2;
                    end
                end
            end
        end
        
        % free carrier absorption
        if all(photon)
            % acoustic 'absorption'
            P(:,24) = W_ac * sqrt(e_ke + photon(1)) .* (2*e_ke + photon(1));
            % 'emission' (the same process, with nearly elastic energy exchange)
            P(:,25) = P(:,24);
            
            % impurity absorption/emission of momentum
            P(:,26) = W_imp * acoth(sqrt(1+photon(1)./e_ke)) ./ sqrt(e_ke);
            P(:,27) = P(:,26);
        end
    end

    function k = ff_k(k,t0,ff,DC,AC)
        % this function does the free flight k propagation
        % note: t - t0 = ff (free flight duration)
        % k(t) = k0 + A(t-t0) + B/w*[cos(wt0) - cos(wt)]
        k = k + DC.*ff; % DC
        if AC(2) ~= 0
            k = k + AC(1)/AC(2)*(cos(AC(2)*t0) - cos(AC(2)*(t0+ff))); % AC
        end
    end


    function q = ff_q(q,k0,t0,ff,DC,AC)
        % this function does the free flight q propagation
        % r(t) = r_0 + k0/m (t-t0) + A/2m (t-t0)^2 + B/mw[cos(wt0)(t-t0) - 1/w(sin(wt) - sin(wt0))
        q = q + k0/m_star.*ff + DC.*ff.^2/(2*m_star);
        if AC(2) ~= 0
            q = q + AC(1)/(m_star*AC(2))*(cos(AC(2)*t0).*ff - 1/AC(2)*(sin(AC(2)*(t0+ff)) - sin(AC(2)*t0)));
        end
    end

    function k = rotate(k,k_old,k_new,beta,phi)
        % function turns beta/phi spherical coordinate rotations into
        % cartesian rotations
        cosb = cos(beta);
        sinb = sin(beta);
        cosp = cos(phi);
        sinp = sin(phi);
        k_xy_n = sqrt(sum(k(:,1:2).^2,2)); % x-y projection of k vector

        % cos(gamma)
        cosy = k(:,3)./k_old; % kz
        siny = k_xy_n./k_old;
        cosy2 = k(:,2)./k_old; % ky
        cosy3 = k(:,1)./k_old; % kx

        % eta
        cosn = k(:,2)./k_xy_n;
        sinn = k(:,1)./k_xy_n;

        % rotate to k'
        k(:,3) = k_new.*(cosb.*cosy + sinb.*cosp.*siny);
        k(:,2) = k_new.*(cosb.*cosy2 - sinb.*cosp.*cosy.*cosn - sinb.*sinp.*sinn);
        k(:,1) = k_new.*(cosb.*cosy3 - sinb.*cosp.*sinn.*cosy + sinb.*sinp.*cosn);
    end
    %% assign the GA variables and rearrange them a little for coding convenience
    totalStart = tic;
    
    % use the twistah!
    rng('shuffle','simdTwister')
    
    % vars is what the GA sends in
    [param,~,~] = mcVary(param,vars);
    
    % combine separate user inputs into one long vector, pad with zeros
    comp_length = max([cellfun(@length,{param.comp_x,param.comp_y,param.comp_z}) 2]);
    comp_x = padarray(param.comp_x,[0 comp_length-length(param.comp_x)],0,'post');
    comp_y = padarray(param.comp_y,[0 comp_length-length(param.comp_y)],0,'post');
    comp_z = padarray(param.comp_z,[0 comp_length-length(param.comp_z)],0,'post');
    imp_length = max(cellfun(@length,{param.imp_x,param.imp_y,param.imp_z}));
    imp_x = padarray(param.imp_x,[0 imp_length-length(param.imp_x)],0,'post');
    imp_y = padarray(param.imp_y,[0 imp_length-length(param.imp_y)],0,'post');
    imp_z = padarray(param.imp_z,[0 imp_length-length(param.imp_z)],0,'post');
    comp = reshape([comp_x comp_y comp_z],[],3);
    imp = reshape([imp_x imp_y imp_z],[],3);
    
    % switches
    doPBC = [param.doPBCy param.doPBCz]; % reflections - PBCx = 0 is electrode

    % param to atomic units
    % convert to meters, then divide by a0
    L = [param.Lx param.Ly param.Lz]*1E-9/param.a0;
    % multiply length by the number of periods
    L = L .* (comp(2,:) + (comp(2,:) == 0));
    % assumed input is in cm^-3, convert to a0^-3
    ni = param.ni*1E6*param.a0^3;
    sim_time = param.sim_time/param.aut;
    % user given field given in V/cm, divide by auE
    % convert from Hz -> a.u. and also add in 2*pi
    E.DC_in = [param.Ex(1) param.Ey(1) param.Ez(1)]/param.auE;
    E.DC = repmat(E.DC_in,param.Ne,1);
    E.AC = [param.Ex(2:3)'.*[1/param.auE;2*pi*param.aut] ...
            param.Ey(2:3)'.*[1/param.auE;2*pi*param.aut] ...
            param.Ez(2:3)'.*[1/param.auE;2*pi*param.aut]];
    % photon parameters
    photon = param.photon .* [1/param.Eh param.wEv*param.a0^3/param.Eh];
    
    % normalize the gradients given to something more accessible for the code
    [comp,c_norm,c_min,c_layout] = grad_normalize(comp,L);
    [imp,i_norm,i_min,i_layout] = grad_normalize(imp,L);
    
    %% MC variables / precalculating scattering curves (for silicon, from Lundstrom, pg 114)
    kT0 = param.k*param.T0/param.Eh; % equilibrium temperature of crystal
    kTe = param.k*param.Te/param.Eh; % initial temperature of electrons
    m_star = 0.26; % conductivity eff mass
    
    % scattering parameters
    Z = 1; % charge
    e0 = 1/(4*pi); % atomic units vacuum permittivity
    kappa = 11.68; % dielectric constant
    %EB = (1/2) * (m_star/kappa^2); % scaled binding energy (half an Eh = 13.6 eV)
    EB = 0.043 / param.Eh; % As: 0.054 / P: 0.045 / Sb: 0.043
    
    % determine if ionized scattering or neutral scattering dominates
    m_dos = 1.08; % density of states effective mass
    Nc = 2 * (m_dos*kT0/(2*pi))^(3/2);
    N_star = Nc/2 * exp(-EB/kT0);
    % n_0 is the free carrier density, i.e. the number of ionized dopants
    % (assuming that that number is much greater than the intrinsic). n_0
    % is used instead of n_i, the dopant concentration, because only
    % ionized dopants can exert a coulombic scattering force
    n_0 = -N_star/2 + sqrt(N_star^2/4 + N_star*ni);
    % some cutoff % of ionization - higher because the ionized scattering
    % rate also decreases with temperature, instead of just looking at % of
    % ionized dopants
    if n_0/ni < 0.25
        unionized = 1;
    else
        unionized = 0;
    end
    fprintf('Unionized?: %g',unionized);
    
    % BH
    beta = sqrt(4*pi*n_0/(kappa*kT0)); % inverse Debye screening length
    eps_b = beta^2/(2*m_star); % beta^2/2m
    P_imp_BH = 2^(5/2)*pi*Z^2/(kappa^2*eps_b^2*sqrt(m_star));
    % CW
    b_imp = (3/(4*pi*n_0))^(1/3); % maximum impact parameter
    e_b = 1/(2*kappa*b_imp);
    P_imp_CW = pi*n_0*Z^2*b_imp^2*sqrt(2/m_star);
    if b_imp > 1/beta
        % if the screening length is smaller than the mean distance between
        % impurities, use BH scattering
        BH = 1;
    else
        % else, use CW scattering
        BH = 0;
    end
    fprintf('\nBH?: %g\n',BH);
    
    % hydrogenic excitation / ionization rates
    a_0 = 1/(2*kappa*e0*EB);
    % inelastic (1 -> 2)
    n = 1; % start level
    n_p = 2; % end level
    x = 1 - (n/n_p)^2;
    if n == 1
        rn = 0.45;
        g0 = 1.1330;
        g1 = -0.4059;
        g2 = 0.07014;
        bn = -0.603;
    elseif n == 2
        rn = 1.94*n^(-1.57);
        g0 = 1.0785;
        g1 = -0.2319;
        g2 = 0.02947;
        bn = 4.0/n - 18.63/n^2 + 36.42/n^3 - 28.09/n^4;
    else
        rn = 1.94*n^(-1.57);
        g0 = 0.9935 + 0.2328/n - 0.1296/n^2;
        g1 = -(0.6282/n - 0.5598/n^2 + 0.5299/n^3);
        g2 = 0.3887/n^2 - 1.181/n^3 + 1.470/n^4;
        bn = 4.0/n - 18.63/n^2 + 36.42/n^3 - 28.09/n^4;
    end
    fnn = 32/(3*sqrt(3)*pi) * n/n_p^3 * x^-3 * (g0 + g1/x + g2/x^2);
    Enn = EB * (1/n^2 - 1/n_p^2);
    rnn = rn * x;
    Ann = 2*n^2*fnn/x;
    Bnn = 4*n^4/n_p^3 * x^-2 * (1 + 4/3*x^-1 + bn*x^-2);
    % inelastic (ionization)
    An = 32/(3*sqrt(3)*pi) * n * sum([g0 g1 g2]./((0:2)+3));
    Bn = 2/3 * n^2 * (5 + bn);

    % acoustic phonon
    E1 = 9.5/param.Eh; % acoustic deformation potential
    ul = 9.04*1E3/param.a0*param.aut; % sound velocity (longitudinal)
    e_u = m_star*ul^2/2; % energy of e- at sound velocity
    rho = 2329*param.a0^3/param.erm; % crystal density
    C_ac = 4*sqrt(e_u)/kT0; % prefactor for x in F(x) and G(x) integrals
    Pe_ac = sqrt(m_star)*kT0^3*E1^2/(sqrt(2)^5*pi*ul^4*rho); % prefactor
    Pkq_ac = m_star*E1^2/(4*pi*rho*ul); % prefactor

    % intervalley scattering (acoustic + optical)
    Zf = [1 1 1 4 4 4]; % number of possible final equivalent valleys
    % g g g f f f (g = across, f = 4 adjacent)
    DtK = [0.5 0.8 11.0 0.3 2.0 2.0]*1E10/param.Eh*param.a0; % coupling constants
    wop = [0.012 0.019 0.062 0.019 0.047 0.059]/param.Eh; % energy of transition
    Nop = 1./(exp(wop/kT0)-1); % optical boltzmann prefactor
    P_op = Zf.*DtK.^2*sqrt(m_star)^3./(sqrt(2)*pi*rho*wop);

    % the energy range must include e_u to capture the change in the integrals
    samples = 400;
    e_range_top = 3/param.Eh; % in Ha
    e_range = [logspace(log10(sqrt(eps)),log10(e_u),samples) logspace(log10(e_u + sqrt(eps)),log10(e_range_top),2*samples)]';
    
    % precalculate acoustic scattering probabilities (as function of
    % electron kinetic energy) for reference, Si's longitudinal e_u is
    % 0.0604 meV, compared to room temperature 25 meV.
    % first column is absorption, second is emission
    Pac = zeros(numel(e_range),2);
    % prefactors
    x1a = C_ac*(sqrt(e_u) - sqrt(e_range(1:samples))); % e < e_u
    x2a = C_ac*(sqrt(e_u) + sqrt(e_range)); % valid for all e
    x2e = C_ac*(sqrt(e_range(samples+1:end))-sqrt(e_u)); % e > e_u
    % calculate the probabilities
    for a = 1:numel(e_range)
        if a <= samples
            Pac(a,1) = Pe_ac/sqrt(e_range(a))*integral(@F1,x1a(a),x2a(a));
        else
            Pac(a,1) = Pe_ac/sqrt(e_range(a))*integral(@F1,0,x2a(a));
            Pac(a,2) = Pe_ac/sqrt(e_range(a))*integral(@G1,0,x2e(a-samples));
        end
    end
    P_ac_abs = griddedInterpolant(e_range,Pac(:,1));
    P_ac_em = griddedInterpolant(e_range,Pac(:,2));
    
    % photon absorption
    W_ac = sqrt(2*m_star)/(3*pi*kappa*e0) * photon(2)/photon(1)^4 * E1^2*kT0/(2*ul^2*rho);
    W_ac_kq = 1/(12*pi*kappa*e0*m_star) * photon(2)/photon(1)^4 * E1^2*kT0/(2*ul^2*rho);
    W_op = sqrt(2*m_star)/(3*pi*kappa*e0) * photon(2)/photon(1)^4 * DtK.^2./(2*wop*rho) .* Zf;
    W_op_kq = 1/(12*pi*kappa*e0*m_star) * photon(2)/photon(1)^4 * DtK.^2./(2*wop*rho) .* Zf;
    W_imp = 1/(3*pi*kappa*e0*sqrt(2*m_star)^3) * photon(2)/photon(1)^4 * Z^2*n_0/(kappa*e0)^2;
    W_imp_kq = 1/(12*pi*kappa*e0*m_star) * photon(2)/photon(1)^4 * Z^2*n_0/(kappa*e0)^2;
    
    % strain - the user supplies two values given in comp_range. the strain
    % effects of the two end points are calculated using linear
    % interpolation formula. then the comp_{x,y,z} functions specify how
    % the strain varies between the two endpoints (d)elta(4)(b)ulk
    d4b = 0.00434*diff(param.comp_range); % get deltaE from absolute energy of d4/d2 fit
    d2b = 0.0286*diff(param.comp_range);
    
    % estimate the highest probability of scattering
    % assume maximal strain for greatest energy difference, 4-fold valley 
    % n_0 is number of ionized dopants
    P_total = scatProb(-0.0246*diff(param.comp_range)*ones(size(e_range)),n_0,e_range,true(size(e_range)));
    P_total = sum(P_total,2);
    P_tot = griddedInterpolant(e_range,P_total);
    
    if param.verbose
        fprintf('Finished parameter initialization & acoustic curve integration: ')
        toc(totalStart)
    end

    %% Propagation and observable collection
    % preallocation & initialization - the electron coordinate and k data
    % structure holds the STARTING values of each free flight. All the
    % resulting data, up until collision, can be found by analytic newton
    % propagation    
    e = struct;
    % electron index goes down the column
    e.q = bsxfun(@times,rand(param.Ne,3),L);
    e.k = normrnd(0,sqrt(m_star*kTe),param.Ne,3);
    e.t = zeros(param.Ne,1);
    e.v = ones(param.Ne,1);
         % electron valley tracker: -1 = 4, 1 = 2 (-1 multiplies the band
         % splitting C so the e- goes downhill in energy, as Efi = (-))
    e.v(rand(param.Ne,1) < exp(-0.0246*diff(param.comp_range)/kTe)) = -1; % some electrons start in higher valleys
    % boltzmann distribute electron's around the device
    for a = 1:3 % x y z
        if comp(2,a) > 0
            bd_q = linspace(0,comp(2,a),1000);
        else
            bd_q = linspace(0,L(a),1000);
        end
        bd = exp(pgrad(comp(:,a),bd_q)*0.0286*diff(param.comp_range)/kTe);
        bd = cumsum(bd/sum(bd));
        for b = 1:param.Ne
            [~,ind_q] = min(abs(bd - rand(1)));
            e.q(b,a) = bd_q(ind_q);
        end
        % handle the 4-fold valley
        bd = exp(pgrad(comp(:,a),bd_q)*0.00434*diff(param.comp_range)/kTe);
        bd = cumsum(bd/sum(bd));
        for b = find(e.v==-1)'
            [~,ind_q] = min(abs(bd - rand(1)));
            e.q(b,a) = bd_q(ind_q);
        end
        
        % if explicitly looking at multiple periods, randomly translate e-
        % by whole periods
        if comp(2,a) > 0 && L(a) > comp(2,a)
            e.q(:,a) = e.q(:,a) + randi(floor(L(a)/comp(2,a)),param.Ne,1)*comp(2,a);
        end
    end
    
    % initialize the observables - time goes down the column
    MQ         = zeros(1,0);
    MQ2        = zeros(1,0);
    MK         = zeros(1,0);
    MK2        = zeros(1,0);
    flux_right = zeros(1,0);
    flux_left  = zeros(1,0);
    en         = zeros(1,0);
    two_occ    = zeros(1,0);
    q_ts       = zeros(1,0); % q_[time series]
    k_ts       = zeros(1,0);
    if param.doX
        % <x>
        MQ = zeros(param.sim_slice+1,3);
        MQ(1,:) = mean(e.q);
    end
    if param.doX2
        MQ2 = zeros(param.sim_slice+1,3);
        MQ2(1,:) = mean(e.q.^2);
    end
    if param.doMV
        % <v>
        MK = zeros(param.sim_slice+1,3);
        MK(1,:) = mean(e.k);
    end
    if param.doMV2
        MK2 = zeros(param.sim_slice+1,3);
        MK2(1,:) = mean(e.k.^2);
    end
    if param.doFlux
        % flux, counted as how many times e- crosses the x boundaries
        flux_right = 0;
        flux_left = 0;
    end
    if param.doEnergy
        en = zeros(param.sim_slice+1,1);
        en(1) = mean(sum(e.k.^2,2));
    end
    if param.doValley
        two_occ = zeros(param.sim_slice+1,1);
        two_occ(1) = sum(e.v==1)/param.Ne;
    end
    if param.doMovie
        % make a movie of dimensions where there is asymmetry or NOT pbc
        movie_ind = 1;
        q_ts = zeros(param.Ne,param.sim_slice,sum(comp(1,:) | imp(1,:)));
        k_ts = zeros(param.Ne,param.sim_slice,4);
        ff_ts = zeros(param.Ne,param.sim_slice);
    end
    % division counter to count how many electrons pass the observation pt
    obs_e = zeros(param.sim_slice+1,1);
    obs_e(1) = 1; % t=0
    
    % estimate gamma
    if BH
        gamma = max(P_tot(e_range));
    else
        gamma = P_tot(1.5*max(sum(e.k.^2,2)/(2*m_star)));
    end
    tau_not = 1/gamma; %t_0
    
    % every amount of this time collect observables
    obs_dt = sim_time/param.sim_slice;

    if param.verbose
        fprintf('Initial tau_not: %.2f\n',tau_not)
        % loop counter
        lc = 0;
        % boundary scattering events
        scatters = zeros(1,2);
        % self scatters
        self_s = 0;
        % photon scatters
        photon_s = 0;
        % maximum energy in simulation
        e_ke_sim_max = max(sum(e.k.^2,2)/(2*m_star));
        % lasso'd e-
        e_lasso = 0;
    end
    
    % random number bank initialize
    rn_size = 1E6;
    rn_account = rand(rn_size,1);
    rn_counter = 0;
    
    % actively propagating electrons (decreases when an e- finishes)
    active_e = param.Ne;
    
    % main while loop
    while active_e > 0 % while there are still electrons propagating
        if param.verbose
            % loop counter
            lc = lc + 1;
        end
        % free flight
        if rn_counter + active_e > rn_size
            rn_account = rand(rn_size,1);
            rn_counter = 0;
        end
        % generate free flight duration
        ff = -tau_not*log(rn_account(rn_counter+1:rn_counter+active_e));
        rn_counter = rn_counter + active_e;
        
        % check if current time + ff passes a marker for observable
        % collection. to determine if the ff crosses a marker, check
        % mod(t) with mod(t+ff) since mod(t+ff) will always roll over
        % to a smaller number (except if time == slice_time...)
        oind = find(mod(e.t,obs_dt) + ff > obs_dt);
        % find the index at which to store observable
        obs_ind = floor((e.t(oind)+ff(oind))/obs_dt);
        % get rid of any observation points over the sim_slice input
        oind(obs_ind > param.sim_slice) = [];
        obs_ind(obs_ind > param.sim_slice) = [];
        % add 1 to offset as the first index stores t=0 conditions
        obs_ind = obs_ind + 1;
        % calculate the ff time to observation point
        obs_ff = ff(oind) - mod(e.t(oind)+ff(oind),obs_dt);
        % split ff into bits
        ff = ff/param.dt_slice;
        % determine at which slice to collect the observable
        obs_ff_ind = floor(obs_ff./ff(oind));
        obs_ff = obs_ff - ff(oind).*obs_ff_ind;

        % loop that propagates sliced up ff
        for a = 0:param.dt_slice-1
            % update the DC fields, only if doing sawtooth - i.e. linear
            for b = 1:3 % x, y, z
                if comp(1,b) == -4
                    if comp(2,b) == L(b) % skip the mod operation
                        cind = e.q(:,b);
                    else
                        cind = mod(e.q(:,b),comp(2,b));
                    end
                    % electrons in d2 band
                    E.DC(cind < comp(3,b)*comp(2,b),b) = d2b/(comp(3,b)*comp(2,b));
                    E.DC(cind > comp(3,b)*comp(2,b),b) = -d2b/((1-comp(3,b))*comp(2,b));
                    % electrons in d4 band
                    E.DC(e.v == -1,b) = E.DC(e.v == -1,b)*d4b/d2b;
                    % add the user supplied DC offset
                    if E.DC_in(b)
                        E.DC(:,b) = E.DC(:,b) + E.DC_in(b);
                    end
                end
            end
            
            collected = 0;
            if any(oind) % direct index
                oind2 = obs_ff_ind == a; % logical index
                if any(oind2) % if we are at a correct slice ..
                    collected = 1;
                    oind3 = oind(oind2);
                    % propagate up to the observation time
                    for b = 1:3
                        e.q(oind3,b) = ff_q(e.q(oind3,b),e.k(oind3,b),e.t(oind3),obs_ff(oind2),...
                                             E.DC(oind3,b),E.AC(:,b));
                        e.k(oind3,b) = ff_k(e.k(oind3,b),e.t(oind3),obs_ff(oind2),...
                                             E.DC(oind3,b),E.AC(:,b));
                    end
                    % collect observables
                    oind4 = obs_ind(oind2);
                    for b = 1:sum(oind2)
                        if param.doX
                            MQ(oind4(b),:) = MQ(oind4(b),:) + e.q(oind3(b),:);
                        end
                        if param.doX2
                            MQ2(oind4(b),:) = MQ2(oind4(b),:) + e.q(oind3(b),:).^2;
                        end
                        if param.doMV % divide by electron mass later
                            MK(oind4(b),:) = MK(oind4(b),:) + e.k(oind3(b),:);
                        end
                        if param.doMV2 || doEnergy
                            % precalculate the sum since MV2 and en use it
                            SK = e.k(oind3(b),:).^2;
                            if param.doMV2 % divide by m^2 later
                                MK2(oind4(b),:) = MK2(oind4(b),:) + SK;
                            end
                            if param.doEnergy % divide by 2*m_star later
                                en(oind4(b)) = en(oind4(b)) + sum(SK);
                            end
                        end
                        if param.doValley
                            two_occ(oind4(b)) = two_occ(oind4(b)) + (e.v(oind3(b))==1);
                        end
                        obs_e(oind4(b)) = obs_e(oind4(b)) + 1;
                    end
                    % increase electron time by free flight amount to obs time
                    e.t(oind3) = e.t(oind3) + obs_ff(oind2);
                    % update the free flight to account for the rest of the
                    % free flight (before it was interrupted for observation)
                    ff(oind3) = ff(oind3) - obs_ff(oind2);
                end
            end
            
            % update q and k with free flights
            e.q(:,1) = ff_q(e.q(:,1),e.k(:,1),e.t,ff,E.DC(:,1),E.AC(:,1));
            e.k(:,1) = ff_k(e.k(:,1),e.t,ff,E.DC(:,1),E.AC(:,1));
            e.q(:,2) = ff_q(e.q(:,2),e.k(:,2),e.t,ff,E.DC(:,2),E.AC(:,2));
            e.k(:,2) = ff_k(e.k(:,2),e.t,ff,E.DC(:,2),E.AC(:,2));
            e.q(:,3) = ff_q(e.q(:,3),e.k(:,3),e.t,ff,E.DC(:,3),E.AC(:,3));
            e.k(:,3) = ff_k(e.k(:,3),e.t,ff,E.DC(:,3),E.AC(:,3));
            % update the electron's time counter
            e.t = e.t + ff;
            if collected
                % revert ff as we are taking constant steps
                ff(oind3) = ff(oind3) + obs_ff(oind2);
            end
            
            % boundary conditions y/z
            reflected = 0;
            if any(doPBC == 0) % any reflections
                old_time = e.t;
            end
            for b = 2:3
                lZ = e.q(:,b) < 0;
                gL = e.q(:,b) > L(b);
                ind = lZ | gL;
                if any(ind)
                    if doPBC(b-1)
                        e.q(ind,b) = mod(e.q(ind,b),L(b));
                    else
                        reflected = 1;
                        % Assume diffuse reflection -  go back in time to 
                        % when the e- hits the boundary, and flip the k 
                        % vector. the next propagation will then start
                        % from the reflected state
                        tr = faster_fzero(e.q(ind,b),e.k(ind,b),e.t(ind),E.DC(ind,b),E.AC(:,b),...
                            L(b)*(e.q(ind,b)>L(b)),-1.1*ff(ind),zeros(sum(ind),1));
                        
                        % reset back to the boundary value (nb: tr is negative)
                        for c = 1:3
                            e.q(ind,c) = ff_q(e.q(ind,c),e.k(ind,c),e.t(ind),tr,E.DC(ind,c),E.AC(:,c));
                            e.k(ind,c) = ff_k(e.k(ind,c),e.t(ind),tr,E.DC(ind,c),E.AC(:,c));
                        end
                        
                        % move time back
                        e.t(ind) = e.t(ind) + tr;

                        % if resetting y-axis, check that the z-coordinate
                        % is also within boundaries, only rotate those
                        if doPBC(2) == 0 && b == 2
                            lZ = e.q(:,3) > 0 & e.q(:,3) < L(3) & lZ;
                            gL = e.q(:,3) > 0 & e.q(:,3) < L(3) & gL;
                            ind = lZ | gL;
                        end
                        
                        if param.specularity > 0
                            % specular reflection
                            spec = rand(sum(ind),1) < param.specularity;
                            ind_spec = ind;
                            ind_spec(ind_spec) = spec;
                            e.k(ind_spec,b) = -e.k(ind_spec,b);
                            % rest scatter diffusely
                            ind(ind) = ~spec;
                            gL = gL & ind;
                        end
                        reflect_knorm = sqrt(sum(e.k(ind,:).^2,2));
                        % diffuse reflection - cosine law
                        theta = asin(sqrt(rand(sum(ind),1)));
                        phi = 2*pi*rand(sum(ind),1);
                        % for e- < 0, use positive hemisphere
                        e.k(ind,b) = reflect_knorm.*cos(theta);
                        e.k(ind,setdiff([1 2 3],b)) = [reflect_knorm.*sin(theta).*cos(phi) reflect_knorm.*sin(theta).*sin(phi)];
                        % e- > L reflect in - direction
                        e.k(gL,b) = -e.k(gL,b);
                        
                        if param.verbose
                            % number of scattering events
                            scatters(b-1) = scatters(b-1) + sum(ind);
                        end
                    end
                end
            end
            % if there was a reflection, finish the propagation
            if reflected
                ff_r = old_time - e.t;
                ind = ff_r > 0;
                for b = 1:3
                    e.q(ind,b) = ff_q(e.q(ind,b),e.k(ind,b),e.t(ind),ff_r(ind),E.DC(ind,b),E.AC(:,b));
                    e.k(ind,b) = ff_k(e.k(ind,b),e.t(ind),ff_r(ind),E.DC(ind,b),E.AC(:,b));
                end
                e.t = old_time;
            end
            % boundary conditions x
            if param.doFlux
                flux_left = flux_left + sum(e.q(:,1) < 0);
                flux_right = flux_right + sum(e.q(:,1) > L(1));
            end
            if param.doPBCx
                % if periodic boundary conditions, an electron that leaves
                % one boundary is modulused back to the other side
                ind = e.q(:,1) < 0 | e.q(:,1) > L(1);
                e.q(ind,1) = mod(e.q(ind,1),L(1));
            else
                % electrodes - if an electron leaves, replace it from the
                % other side with a thermal electron
                ind = e.q(:,1) > L(1);
                e.q(ind,1) = 0;
                e.q(ind,2) = rand(sum(ind),1)*L(2);
                e.q(ind,3) = rand(sum(ind),1)*L(3);
                e.k(ind,:) = normrnd(0,sqrt(m_star*kTe),sum(ind),3);
                e.k(ind,1) = abs(e.k(ind,1));
                % since the old e- left, the new one is frozen until a new ff is calculated
                ff(ind) = 0;
                ind = e.q(:,1) < 0;
                e.q(ind,1) = L(1);
                e.q(ind,2) = rand(sum(ind),1)*L(2);
                e.q(ind,3) = rand(sum(ind),1)*L(3);
                e.k(ind,:) = normrnd(0,sqrt(m_star*kTe),sum(ind),3);
                e.k(ind,1) = -abs(e.k(ind,1));
                ff(ind) = 0;
            end
        end

        % compositional gradient
        if any(comp(1,:)~=0)
            C = pgrad(comp(:,1),e.q(:,1)) + pgrad(comp(:,2),e.q(:,2)) + pgrad(comp(:,3),e.q(:,3));
            % normalize
            C = (C - c_min)*c_norm;
            % set true user given range 0 -> 1
            C = C * diff(param.comp_range) + param.comp_range(1);
            C = 0.0246*C; % 0.67x (eV)(% Ge)
            C = e.v.*C; % depending on valley, -1 or 1
        else
            C = 0.0246*param.co*e.v;
        end
        % impurity gradient
        if any(imp(1,:)~=0)
            % get the gradient profile based on the e- location
            I = pgrad(imp(:,1),e.q(:,1)) + pgrad(imp(:,2),e.q(:,2)) + pgrad(imp(:,3),e.q(:,3));
            % normalize gradient to 0 -> 1 (or user max, if rectangles)
            I = (I - i_min)*i_norm;
            % normalize to be within user supplied range min -> max
            I = I * diff(param.imp_range) + param.imp_range(1);
            if unionized
                % neutral scattering - no screening
                I = I * (ni - n_0);
            else
                error(['Need to calculate the true electron concentration as a function of coord, '...
                'as it affects the scattering rate via screening (i.e. Poisson self consistent solver)'])
            end
        else
            % uniform scattering
            if unionized
                I = ni - n_0;
            else
                I = n_0;
            end            
        end

        % electron kinetc energy hbar^2k^2/2m* (parabolic band)
        e_ke = sum(e.k.^2,2)/(2*m_star);
        e_knorm = sqrt(2*m_star*e_ke);
        % use random numbers to determine the mechanism of scattering:
        % compare to scattering probabilities (by cum. sum), which is
        % stored as 0 (implied) -> P1 -> P1+P2 -> ... -> P1+...+Pn -> Gamma
        P = cumsum(scatProb(C,I,e_ke,e.v==-1),2);
        % random numbers from 0 -> gamma
        if rn_counter + active_e > rn_size
            rn_account = rand(rn_size,1);
            rn_counter = 0;
        end
        % vectorized find - max returns first index (1, since binary)
        % alg_filt is 0 or 1, and 0 cancels rows with self scattering (all
        % zeros, so max thinks 0 is max value and returns index = 1)
        [alg_filt,scat_choice] = max(bsxfun(@lt,rn_account(rn_counter+1:rn_counter+active_e)*gamma,P),[],2);
        scat_choice = alg_filt.*scat_choice;
        rn_counter = rn_counter + active_e;
        if param.verbose
            % amount of self scatters
            self_s = self_s + sum(scat_choice==0);
        end
        % new e- from electrodes self scatter, i.e. do not scatter
        if param.doPBCx == 0
            scat_choice(ff==0) = 0;
        end
        
        % collect sim_slice number of frames
        if param.doMovie && any(e.t > param.doMovie*sim_time) && movie_ind <= param.sim_slice
            for a = 1:size(q_ts,3)
                q_ts(:,movie_ind,a) = e.q(:,a);
            end
            k_ts(:,movie_ind,1:3) = e.k;
            % scattering information
            k_ts(:,movie_ind,4) = scat_choice;
            ff_ts(:,movie_ind) = e.t;
            movie_ind = movie_ind + 1;
        end

        % generate new state
        for a = 1:45
            % find index of all e- that have scattered with mechanism 'a'
            % the scattering choice of 0 corresponds to self scattering
            % - nothing happens here
            ind = scat_choice == a;
            if any(ind) % if any electron has scattered with mechanism a
                scat_knorm = e_knorm(ind);
                n_scat = sum(ind);
                explicit_q = 0;
                switch a
                    case 1
                        % impurity scattering
                        if unionized
                            % isotropic scattering
                            beta = acos(1-2*rand(n_scat,1));
                        else
                            rns = rand(n_scat,1);
                            if BH
                                beta = acos(1 - 2*(1-rns)./(1+4*rns.*e_ke(ind)/eps_b));
                            else
                                beta = acos(((e_ke(ind)/e_b).^2.*rns - 1) ./ ((e_ke(ind)/e_b).^2.*rns + 1));
                            end
                        end
                        % energy isn't changed, so |k| (i.e. scat_knorm) doesn't change
                    case 2
                        % inelastic neutral impurity excitation, 1 -> 2
                        % rotation angle
                        beta = acos(1-2*rand(n_scat,1));
                        
                        % change energy by excitation energy
                        scat_knorm = sqrt(scat_knorm.^2 - 2*m_star*Enn);
                    case 3
                        % inelastic neutral impurity ionization
                        beta = acos(1-2*rand(n_scat,1));
                        
                        % change energy by ionization energy
                        scat_knorm = sqrt(scat_knorm.^2 - 2*m_star*EB/n^2);
                    case 4
                        explicit_q = 1;
                        % acoustic absorption
                        q = zeros(n_scat,1); % magnitude of q vector

                        % maximum q value (absorption + backwards scattering)
                        q_max = 2*(scat_knorm + m_star*ul);
                        % minimum value of |q| cannot be less than 0
                        q_min = max(2*(-scat_knorm + m_star*ul),0);

                        % f(random) will always be < ac_max, as it is the max value here
                        f_max = Pkq_ac./scat_knorm .* q_max.^2./(exp(q_max*ul/(kT0))-1);

                        % choose a value of q by the rejection technique. pick a random
                        % q_r, select if r*max(f in region) < f(q_r)
                        for b = 1:n_scat
                            duple_scaler = [q_max(b)-q_min(b);f_max(b)];
                            duple = rand(2,1).*duple_scaler + [q_min(b);0]; % q_min <-> q_max, and rC
                            while Pkq_ac/scat_knorm(b) * duple(1)^2/(exp(duple(1)*ul/kT0)-1) < duple(2)
                                % keep going while rC > f(x_r)
                                duple = rand(2,1).*duple_scaler + [q_min(b);0];
                            end
                            q(b) = duple(1); % save the q value
                        end

                        % obtain angle between q and k (error in Jacobini&Reggiani)
                        beta = acos((-q/2 + m_star*ul)./scat_knorm);

                        % before = after, k + q = k'
                        % should confirm that the vector k' norm from subtraction
                        % is the same as that found by energy conservation (it is)
                        % sqrt(e_knorm.^2 + 2*m_star*q*ul) == sqrt(scat_kx.^2 + scat_ky.^2)
                    case 5
                        explicit_q = -1;
                        % acoustic emission
                        q = zeros(n_scat,1);

                        % maximum q value: emission takes all k, minus phonon quanta
                        % minimum q value is 0, as -2*(e_knorm + m_star*ul)
                        % can never be greater than 0  
                        q_max = 2*(scat_knorm - m_star*ul);
                        f_max = Pkq_ac./scat_knorm .* q_max.^2.*(1./(exp(q_max*ul/(kT0))-1)+1);

                        for b = 1:n_scat
                            duple_scaler = [q_max(b);f_max(b)];
                            duple = rand(2,1).*duple_scaler;
                            while Pkq_ac/scat_knorm(b) * duple(1)^2 * (1/(exp(duple(1)*ul/(kT0))-1)+1) < duple(2)
                                duple = rand(2,1).*duple_scaler;
                            end
                            q(b) = duple(1); % save the q value
                        end

                        % obtain angle between q and k
                        beta = acos((q/2 + m_star*ul)./scat_knorm);

                        % before = after, k - q = k'
                        % sqrt(e_knorm.^2 - 2*m_star*q*ul) == sqrt(scat_kx.^2 + scat_ky.^2)
                    case {6,7,8}
                        % optical absorption - g-type scattering
                        beta = acos(1-2*rand(n_scat,1));
                        scat_knorm = sqrt(scat_knorm.^2 + 2*m_star*(wop(a-5)));
                    case {9,10,11}
                        % optical emission - g type scattering
                        beta = acos(1-2*rand(n_scat,1));
                        scat_knorm = sqrt(scat_knorm.^2 - 2*m_star*(wop(a-8)));
                    case {12,13,14}
                        % optical absorption - f-type scattering
                        beta = acos(1-2*rand(n_scat,1));
                        scat_knorm = sqrt(scat_knorm.^2 + 2*m_star*(wop(a-8)-C(ind)));
                        % electron scatters to a different valley
                        e.v(ind) = -e.v(ind);
                    case {15,16,17}
                        % optical emission - f type scattering
                        beta = acos(1-2*rand(n_scat,1));
                        scat_knorm = sqrt(scat_knorm.^2 - 2*m_star*(wop(a-11)+C(ind)));
                        % electron scatters to a different valley
                        e.v(ind) = -e.v(ind);
                    case {18,19,20}
                        % optical absorption - f-type 4->4
                        beta = acos(1-2*rand(n_scat,1));
                        scat_knorm = sqrt(scat_knorm.^2 + 2*m_star*(wop(a-14)));
                    case {21,22,23}
                        % optical emission - f-type 4->4
                        beta = acos(1-2*rand(n_scat,1));
                        scat_knorm = sqrt(scat_knorm.^2 - 2*m_star*(wop(a-17)));
                    case {24,25}
                        % photon + acoustic (elastic approx)
                        q = zeros(n_scat,1);
                        q_min = scat_knorm.*(sqrt(1+photon(1)./e_ke(ind)) - 1);
                        q_max = scat_knorm.*(sqrt(1+photon(1)./e_ke(ind)) + 1);
                        f_max = W_ac_kq*q_max.^3./scat_knorm;
                        for b = 1:n_scat
                            duple_scaler = [q_max(b)-q_min(b);f_max(b)];
                            duple = rand(2,1).*duple_scaler + [q_min(b);0]; % q_min <-> q_max, and rC
                            while W_ac_kq*duple(1)^3/scat_knorm(b) < duple(2)
                                duple = rand(2,1).*duple_scaler + [q_min(b);0];
                            end
                            q(b) = duple(1);
                        end
                        if a == 24 % "add" momentum
                            explicit_q = 1;
                            beta = acos((-q/2 + m_star*photon(1)./q)./scat_knorm);
                        else
                            explicit_q = -1; % "substract" momentum
                            beta = acos((q/2 - m_star*photon(1)./q)./scat_knorm);
                        end
                        photon_s = photon_s + n_scat;
                    case {26,27}
                        % photon + impurity                 
                        q = zeros(n_scat,1);
                        q_min = scat_knorm.*(sqrt(1+photon(1)./e_ke(ind)) - 1);
                        q_max = scat_knorm.*(sqrt(1+photon(1)./e_ke(ind)) + 1);
                        f_max = W_imp_kq./(scat_knorm.*q_max);
                        for b = 1:n_scat
                            duple_scaler = [q_max(b)-q_min(b);f_max(b)];
                            duple = rand(2,1).*duple_scaler + [q_min(b);0]; % q_min <-> q_max, and rC
                            while W_imp_kq/(scat_knorm(b)*duple(1)) < duple(2)
                                duple = rand(2,1).*duple_scaler + [q_min(b);0];
                            end
                            q(b) = duple(1);
                        end
                        if a == 26 % "add" momentum
                            explicit_q = 1;
                            beta = acos((-q/2 + m_star*photon(1)./q)./scat_knorm);
                        else
                            explicit_q = -1; % "substract" momentum
                            beta = acos((q/2 - m_star*photon(1)./q)./scat_knorm);
                        end
                        photon_s = photon_s + n_scat;
                    case {28,29,30}
                        % optical absorption - g-type scattering
                        explicit_q = 1;
                        W_op_ind = a-27;
                        omega_pm = photon(1) + wop(W_op_ind);
                    case {31,32,33}
                        % optical emission - g type scattering
                        explicit_q = -1;
                        W_op_ind = a-30;
                        omega_pm = photon(1) - wop(W_op_ind);
                    case {34,35,36}
                        % optical absorption - f-type scattering
                        explicit_q = 1;
                        W_op_ind = a-30;
                        omega_pm = photon(1) + wop(W_op_ind) - C(ind);
                        % electron scatterings to a different valley
                        e.v(ind) = -e.v(ind);
                    case {37,38,39}
                        % optical emission - f type scattering
                        explicit_q = -1;
                        W_op_ind = a-33;
                        omega_pm = photon(1) - wop(W_op_ind) - C(ind);
                        % electron scatters to a different valley
                        e.v(ind) = -e.v(ind);
                    case {40,41,42}
                        % optical absorption - f-type 4->4
                        explicit_q = 1;
                        W_op_ind = a-36;
                        omega_pm = photon(1) + wop(W_op_ind);
                    case {43,44,45}
                        % optical emission - f-type 4->4
                        explicit_q = -1;
                        W_op_ind = a-39;
                        omega_pm = photon(1) - wop(W_op_ind);
                end
                
                if a>=28 && a<=45
                    % photon + intervalley ("optical") phonon
                    q = zeros(n_scat,1);
                    q_min = scat_knorm.*(sqrt(1+omega_pm./e_ke(ind)) - 1);
                    q_max = scat_knorm.*(sqrt(1+omega_pm./e_ke(ind)) + 1);
                    f_max = W_op_kq(W_op_ind)*q_max.^3./scat_knorm*(Nop(W_op_ind) + 1/2 - explicit_q/2);
                    for b = 1:n_scat
                        duple_scaler = [q_max(b)-q_min(b);f_max(b)];
                        duple = rand(2,1).*duple_scaler + [q_min(b);0]; % q_min <-> q_max, and rC
                        while W_op_kq(W_op_ind)*duple(1)^3/scat_knorm(b)*(Nop(W_op_ind) + 1/2 - explicit_q/2) < duple(2)
                            duple = rand(2,1).*duple_scaler + [q_min(b);0];
                        end
                        q(b) = duple(1);
                    end
                    beta = acos((-explicit_q*q/2 + explicit_q*m_star*omega_pm./q)./scat_knorm);
                    photon_s = photon_s + n_scat;
                end
                
                % the azimuth is chosen randomly
                phi = 2*pi*rand(n_scat,1);
                if explicit_q
                    % if q is explicitly generated, rotate q to the new direction
                    % k' = k + q
                    e.k(ind,:) = e.k(ind,:) + explicit_q*rotate(e.k(ind,:),e_knorm(ind),q,beta,phi);
                else
                    % isotropic scattering
                    % rotate from k to k' and renormalize (if energy changes)
                    e.k(ind,:) = rotate(e.k(ind,:),e_knorm(ind),scat_knorm,beta,phi);
                end
            end
        end
        
        % BH scattering has a very strong peak at low energies, and 
        % updating the estimate of gamma is pointless
        e_ke_max = max(sum(e.k.^2,2)/(2*m_star));
        if ~BH
            gamma = P_tot(1.5*e_ke_max);
            tau_not = 1/gamma; %t_0
        end
        if param.verbose
            if e_ke_max > e_ke_sim_max
                e_ke_sim_max = e_ke_max;
            end
        end
        
        % if any electrons are running away .. lasso them!
        if any(e_ke > e_range_top)
            ind = e_ke > e_range_top;
            % reset them back to thermal e-
            e.k(ind,:) = normrnd(0,sqrt(m_star*kTe),sum(ind),3);
            % count them
            e_lasso = e_lasso + sum(ind);
        end
        
        % any electrons that have been simulated to the user given time
        % are deleted here - some time is wasted on scattering, but as long
        % as they are deleted before the next loop's free flight it is OK
        ind = e.t > sim_time;
        if any(ind)
            if active_e == param.Ne
                if param.verbose
                    fprintf('tau_not at first deletion: %.2f\n',tau_not)
                end
            end
            e.q(ind,:) = [];
            e.k(ind,:) = [];
            e.t(ind) = [];
            e.v(ind) = [];
            E.DC(ind,:) = [];
            active_e = active_e - sum(ind);
        end
    end
    
    if param.verbose
        fprintf('Maximum energy reached: %.3f eV\n',e_ke_sim_max*param.Eh)
        fprintf('Number of y/z scatters occured: %u %u\n',scatters)
        fprintf('Number of while loop ran: %u\n',lc)
        fprintf('Number of self scatters: %u\n',self_s)
        fprintf('Number of photon scatters: %u\n',photon_s)
        fprintf('Number of e- lassos: %u\n',e_lasso)
        fprintf('Final tau_not: %.2f\n',tau_not)
    end
    
    %% output
    
    % dividing by obs_e instead of Ne helps ensure that an electron with an
    % abnormally long free flight that goes over two observation periods
    % does not count towards averaging the slice that it missed
    if param.doX
        MQ = bsxfun(@rdivide,MQ,obs_e) * param.a0/1E-9;
    end
    if param.doX2
        MQ2 = bsxfun(@rdivide,MQ2,obs_e) * (param.a0/1E-9)^2;
    end
    if param.doMV
        MK = bsxfun(@rdivide,MK,obs_e)/m_star;
    end
    if param.doMV2
        MK2 = bsxfun(@rdivide,MK2,obs_e)/m_star^2;
    end
    if param.doEnergy
        en = en./obs_e/(2*m_star);
    end
    if param.doValley
        two_occ = two_occ./obs_e;
    end
    
    if param.doAvg
        % if doAvg is enabled, any observable is returned as an average
        if param.doX % <x>
            MQ = mean(MQ);
        end
        if param.doX2 % <x^2>
            MQ2 = mean(MQ2);
        end
        if param.doMV % <hk/m>
            MK = mean(MK);            
        end
        if param.doMV2 % <v^2> = <(hk/m)^2>
            MK2 = mean(MK.^2);
        end
        if param.doEnergy % p^2/2m = (hk)^2/2m
            en = mean(en);
        end
        if param.doValley
            two_occ = mean(two_occ);
        end
        fit_value = [MK MK2 MQ MQ2 flux_right flux_left en two_occ];
        
        if param.doEff
            % efficiency defined as <v>^2/<v^2>; if doing a GA calculation, 
            % return absolute value negative (find most negative value), 
            % override other fit_values
            fit_value = -abs(MK^2/MK2);
            
            % if doing GA, save the output
            if param.saveIntermediate
                temp = [vars MK MK2 MQ MQ2 flux_right flux_left en fit_value two_occ];
                % save intermediate data to a text file
                f = fullfile('defectOutput',['dynOpt_intermediate_' param.job_title '_' param.tor '.int']);
                save(f,'temp','-ascii','-append'); 
            end
        end
    else
        % return full output (x,y,z,...)
        fit_value.MQ         = MQ;
        fit_value.MQ2        = MQ2;
        fit_value.MK         = MK;
        fit_value.MK2        = MK2;
        fit_value.flux_left  = flux_left;
        fit_value.flux_right = flux_right;
        fit_value.en         = en;
        fit_value.obs_e      = obs_e;
        fit_value.two_occ    = two_occ;
        fit_value.q_ts       = q_ts * param.a0/1E-9; % coordinate
        fit_value.k_ts       = k_ts; % raw k-value
        fit_value.ff_ts      = ff_ts;
        if param.verbose
            fit_value.max_en     = e_ke_sim_max;
            fit_value.yz_scatter = scatters;
            fit_value.lc         = lc;
            fit_value.self_s     = self_s;
            fit_value.photon_s   = photon_s;
            fit_value.lassos     = e_lasso;
            fit_value.end_tau    = tau_not;
        end
        if strcmp(getenv('PBS_JOBID'),'') % if running on a local computer
            fit_value.c_layout   = ((c_layout-c_min)*c_norm) * diff(param.comp_range) + param.comp_range(1);
            fit_value.i_layout   = ((i_layout-i_min)*i_norm) * diff(param.imp_range) + param.imp_range(1);
        end
    end
    %% toc
    fprintf('Finished a manchego-carlo fitness run: %.2f seconds\n',toc(totalStart));
end