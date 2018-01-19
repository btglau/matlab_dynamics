function [ Vt ] = propagate_1D_Vt(t,param)
%% PROPAGATE_1D_VT take a time and input structure and return a value from 0
%to 1 that controls the time depedence of the potential. Basically, return
%V(t) part of V(x,t).

    Vt = zeros(param.unique_tiles,length(t));
    i_ = 1; % variable tracking fourier norm
    
    for a = 1:param.unique_tiles
        ms = param.time_matrix(param.time_ind(a):param.time_ind(a+1)-1);
        switch ms(1)
            case -1
                % periodic: mod(time,period)
                per_t = mod(t,sum(ms(2:5)));
                switch ms(6)
                    % phases
                    case 1
                        % start with 0
                        ms_cs = cumsum(ms(2:5));
                        Vt(a,per_t < ms_cs(1)) = 0;
                            temp_index = ms_cs(1) <= per_t & per_t < ms_cs(2);
                        Vt(a,temp_index) = 1/ms(3)*(per_t(temp_index)-ms_cs(1));
                        Vt(a,ms_cs(2) <= per_t & per_t < ms_cs(3)) = 1;
                            temp_index = ms_cs(3) <= per_t & per_t < ms_cs(4);
                        Vt(a,temp_index) = 1 - 1/ms(5)*(per_t(temp_index)-ms_cs(3));
                    case 2
                        % start with rise
                        ms_cs = cumsum([ms(3:5) ms(2)]);
                            temp_index = per_t < ms_cs(1);
                        Vt(a,temp_index) = 1/ms(3)*per_t(temp_index);
                        Vt(a,ms_cs(1) <= per_t & per_t < ms_cs(2)) = 1;
                            temp_index = ms_cs(2) <= per_t & per_t < ms_cs(3);
                        Vt(a,temp_index) = 1 - 1/ms(5)*(per_t(temp_index)-ms_cs(2));
                        Vt(a,ms_cs(3) <= per_t & per_t < ms_cs(4)) = 0;
                    case 3
                        % start with 1
                        ms_cs = cumsum([ms(4:5) ms(2:3)]);
                        Vt(a,per_t < ms_cs(1)) = 1;
                            temp_index = ms_cs(1) <= per_t & per_t < ms_cs(2);
                        Vt(a,temp_index) = 1 - 1/ms(3)*(per_t(temp_index)-ms_cs(1));
                        Vt(a,ms_cs(2) <= per_t & per_t < ms_cs(3)) = 0;
                            temp_index = ms_cs(3) <= per_t & per_t < ms_cs(4);
                        Vt(a,temp_index) = 1/ms(5)*(per_t(temp_index)-ms_cs(3));
                    case 4
                        % start with fall
                        ms_cs = cumsum([ms(5) ms(2:4)]);
                            temp_index = per_t < ms_cs(1);
                        Vt(a,temp_index) = 1 - 1/ms(5)*per_t(temp_index);
                        Vt(a,ms_cs(1) <= per_t & per_t < ms_cs(2)) = 0;
                            temp_index = ms_cs(2) <= per_t & per_t < ms_cs(3);
                        Vt(a,temp_index) = 1/ms(3)*(per_t(temp_index)-ms_cs(2));
                        Vt(a,ms_cs(3) <= per_t & per_t < ms_cs(4)) = 1;
                end
            case -2
                % build time argument
                per_t = pi*t/ms(2);
                
                % chop out the -2 and tau to make indexing in loop easier
                ms = ms(3:end);
                
                % bsxfun if time vector is long
                if length(t) > 1
                    % fourier series by bsxfun
                    Vt(a,:) = sum(bsxfun(@times,ms(1:2:end),sin(bsxfun(@plus,bsxfun(@times,(1:length(ms)/2),per_t),ms(2:2:end))).^2),2);
                else % a single time point (no bsxfun)
                    Vt(a,1) = sum(ms(1:2:end).*sin((1:length(ms)/2)*per_t + ms(2:2:end)).^2);
                end
                
                % more than one fourier term
                if length(ms) > 2
                    % normalize this row of vt so it goes between 0 and 1
                    Vt(a,:) = Vt(a,:)/param.fourier_norm(i_);
                    i_ = i_ + 1;
                end
        end
    end
end

