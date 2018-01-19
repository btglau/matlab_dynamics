function [T,kxs,Tr] = initializeT_1D(param)

% CONVENTION:
% Rows -> Y
% Cols -> X

% Initialize T
% Loop over a NxN matrix with correct k-values. Check dimension.
% 1 -> N/2 + 1;   N/2+2 -> N   (index in for loop in matlab)
% 0 -> N/2    ; -(N/2-1) -> -1 (on paper shifted matrix to line up with FFT)
% conversion: (-1) and -(N+1)
% The conversion is to make the values of k match up with the FFT algorithm
    
    % Create a vector of wavenumbers (shifted to match FFT indexing)
    if mod(param.Nx,2) == 0 % even grid points
        kxs = param.dkx*[0:1:param.Nx/2 (-param.Nx/2+1):1:-1].';
    else % odd number of grid points
        kxs = param.dkx*[0:1:param.Nx/2 -param.Nx/2:1:-1].';
    end
    
    % T = k^2/(2m_eff)
    T = kxs.^2/(2*param.effMass);
    
    Tr = [];
    
    if fix(abs(param.bias)) == 9
        % create the real space representation of T (Tannor, periodic)
        % off diagonal is j-i, where <i|T|j>. The actual index start and
        % end doesn't matter since it's just subtraction (11.172)
       Tr = bsxfun(@minus,(0:param.Nx-1),(0:param.Nx-1).');
       if mod(param.Nx,2) == 0 % even grid points
           % off diagonal elements (2/N^2)*(-1)^j-i / sin^2(pi(j-i)/N)
           Tr = (2/param.Nx^2)*((-1).^Tr)./sin(pi*Tr/param.Nx).^2;
           % handle the diagonal (1+2/N^2)/3
           Tr(1:param.Nx+1:param.Nx^2) = (1 + 2/param.Nx^2)/3;
           % constants (K^2/2m)
           Tr = (pi/param.dx)^2/(2*param.effMass)*Tr;
       else % odd grid points
           Tr = (2/param.Nx^2)*((-1).^Tr).*cos(pi*Tr/param.Nx)./sin(pi*Tr/param.Nx).^2;
           Tr(1:param.Nx+1:param.Nx^2) = (1 + 1/param.Nx^2)/3;
           Tr = (pi/param.dx)^2/(2*param.effMass)*Tr;
       end
    end
end