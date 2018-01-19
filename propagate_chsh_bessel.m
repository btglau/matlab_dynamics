function [chshBessel,N] = propagate_chsh_bessel(alpha)

% This function returns a vector containing the bessel coefficients
% needed to converge the chebyshev expansion - the 'a_n' (the factor of 2
% is built in)

    % Determine the convergence criteria 
    
    if isreal(alpha)
        N = ceil(alpha); % This is the max number of polynomials needed for good convergence

        while 1
            if besselj(N,alpha) < sqrt(eps) % Find the N that corresponds to the Bessel coefficient under some limit
                N = N+1;
                break
            end
            N = N+1;
        end

        %chshBessel = zeros(N,1);
        %chshBessel(1) = besselj(0,alpha);
        %chshBessel(2:N) = 2*besselj(1:N-1,alpha);
        
        chshBessel = [besselj(0,alpha) 2*besselj(1:N-1,alpha)];
    end
    
    % note that bad things happen with imaginary bessel functions -
    % modified bessel functions of the first kind are exponentially growing
    % and *large* imaginary time steps will cause divergence of the expansion
    if ~isreal(alpha) % modified bessel functions of the first kind for imag args
        alpha = -1i*alpha;
        N = ceil(abs(alpha)); % This is the max number of polynomials needed for good convergence

        while 1
            if besseli(N,alpha) < sqrt(eps) % Find the N that corresponds to the Bessel coefficient with 1E-10 (very small contribution)
                N = N+1;
                break
            end
            N = N+1;
        end

        %chshBessel = zeros(N,1);
        %chshBessel(1) = besseli(0,alpha);
        %chshBessel(2:N) = 2*besseli(1:N-1,alpha);
        
        chshBessel = [besseli(0,alpha) 2*besseli(1:N-1,alpha)];
    end
end