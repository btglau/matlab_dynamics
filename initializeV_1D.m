function [V,imagV,param] = initializeV_1D(param)
%% Initialize 1D ratchet potential for one period

    % column vector of potentials
    V = zeros(param.Nx,param.unique_tiles);
    ru = param.space_matrix(1); % repeat units
    
    for a = 1:param.unique_tiles
        ms = param.space_matrix(param.space_ind(a):param.space_ind(a+1)-1);
        switch ms(1)
            case -1
                % length of rectangle = Nx/number of rectangles
                rect_length = floor(param.Nx/length(ms(2:end)));
                
                % correct further if more than 1 repeat unit is specified
                rect_length = floor(rect_length/ru);
                
                % create the repeating ms structure based upon repeats
                ms = repmat(ms(2:end-1),1,ru);
                
                for b = 1:length(ms)
                    V((b-1)*rect_length+1:(b)*rect_length,a) = ms(b);
                end
                
                % pad the leftover zeros as the last rectangle
                V(b*rect_length:end,a) = ms(b);
            case -2
                % 2 pi x/(L/repeat units)
                X = 2*pi * param.X/(param.Lx/ru);
                
                % chop out the -2 to make indexing in loop easier
                ms = ms(2:end);
                
                % fourier series
                V(:,a) = sum(bsxfun(@times,ms(1:2:end),sin(bsxfun(@plus,bsxfun(@times,(1:length(ms)/2),X),ms(2:2:end)))),2);
        end
    end
    
    % ---------------------------------------------------------------------
    % NIP
    
    % null output
    imagV = [];
    
    if param.doNIP

        % preallocate
        imagV = zeros(param.Nx,1);

        % If the first element of absorb is actually a number, generate a
        % linear absorbing potential between absorb(2) and absorb(3) with
        % slope absorb (1)
        if param.absorbx(1) < 777 && param.absorbx(1) > 1
            % Add absorbing barrier, x
            ab = [param.absorbx(1) 0]; % create the linear polynomial

            for a = param.absorbx(2):param.absorbx(3) % do the right barrier, start->stop
                imagV(a) = imagV(a) -1i*polyval(ab,a-absorbx(2));
            end
            for a = param.Nx-param.absorbx(3)+1:param.Nx-param.absorbx(2) % do the left barrier, reverse from stop->start
                imagV(a) = imagV(a) -1i*polyval(ab,param.Nx-param.absorbx(2)-a);
            end
        end

        % Transmissionless 'elliptical NIP'
        c = 2.62;
        rx = param.absorbx(3)-param.absorbx(2);
        x = c*((param.absorbx(2):param.absorbx(3)).'-param.absorbx(2))/rx;

        if param.absorbx(1) == 888 || param.absorbx(1) == 999
            imagV(param.absorbx(2):param.absorbx(3)) = imagV(param.absorbx(2):param.absorbx(3)) - 2i*(pi/rx)^2*(4./(c-x).^2 + 4./(c+x).^2 - 8/(c)^2);
        end
        if param.absorbx(1) == 777 || param.absorbx(1) == 999
            imagV((param.Nx-param.absorbx(2)+1):-1:(param.Nx-param.absorbx(3)+1)) = imagV((param.Nx-param.absorbx(2)+1):-1:(param.Nx-param.absorbx(3)+1)) - 2i*(pi/rx)^2*(4./(c-x).^2 + 4./(c+x).^2 - 8/(c)^2);
        end
    end
end