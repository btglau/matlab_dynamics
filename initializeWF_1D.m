function [Psi] = initializeWF_1D(param,KE)

% This function creates a product gaussian wave packet

    switch param.doBath(1)
        case 1
            % (S)SH, many column spinor wavefunction
            Psi = zeros(param.Nx,2^param.bathN);
        otherwise
            % bath is off, single column wavefunction
            Psi = zeros(param.Nx,1);
    end

% Initialize the UNNORMALIZED wavefunction

    % for a given kinetic energy find the closest possible allowed k-space
    % value in the (+) part of the discrete spectrum
    %{
    kspace = param.dkx*(0:(param.Nx/2));
    [~,ik] = min(abs(kspace-sqrt(abs(2*param.effMass*KE))));
    kx = kspace(ik);
    %}
    
    % ignore above code: any k-space value should be able to be represented
    % as a linear combination of the discrete k-space spectrum (that's what
    % the fourier transform is for)
    kx = sqrt(2*param.effMass*KE);
    
    % change the sign of kx based upon the direction bias selected
    switch param.bias
        case {1,4}
            kx = -kx; % left running wavepacket
        case {2,5}
            % spreading wavepacket; kx stays positive
        case {3,6}
            % kx stays positive for right running wavepacket
    end
    
    switch fix(abs(param.bias)) % look at integer part
        case {1,3}
            % gaussian plane wave packet
            Psi(:,1) = exp(-(param.X-param.x0).^2/(4*param.sigmax^2)).*exp(1i*kx*param.X);
        case 2
            % gaussian cosine, centered at x0
            Psi(:,1) = exp(-(param.X-param.x0).^2/(4*param.sigmax^2)).*cos(kx*(param.X-param.x0));
        case {4,6}
            % plane wave
            Psi(:,1) = exp(1i*kx*(param.X-param.x0));
        case 5
            % left + right plane wave
            Psi(:,1) = cos(kx*(param.X-param.x0));
        case 7 % uniform amplitude, random phase wavefunction
            Psi(:,1) = exp(1i*random('unif',0,1,[1 param.Nx]))/param.Nx;
        case 8 % gaussian amplitude, random phase
            Psi(:,1) = exp(-(param.X-param.x0).^2/(4*param.sigmax^2)).*exp(1i*random('unif',0,1,[1 param.Nx]));
        case 9 % eigenfunction calc
            pot_avg = abs(param.bias) - fix(abs(param.bias));
            [Psi,eval] = eig(param.Tr + diag(param.V*pot_avg),'vector');
            Psi = [Psi eval];
    end
    
    % condition Psi
    Psi(abs(Psi)<sqrt(eps)) = 0;
    
    % normalize the wavefunction
    Psi(:,1) = Psi(:,1)/norm(Psi(:,1),2);
end