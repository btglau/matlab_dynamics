function [nan_var,lb,ub] = split_input(str)
% this function takes the string input, which is a list of values separated
% by 'r', each separated by a ',', and splits it into three vectors, lb
% (left value), ub (right value), and a numeric vector that has NaN
% where a range was given, or the value where a single value was given

%%%%% don't fuck around with special characters (i.e. shift+numeric). Linux
%%%%% will scrub them or do weird things with them.

    % force string to be a row vector
    str = reshape(str,1,[]);

    if ischar(str)
        % remove brackets from the beginning and end (needed for MATLAB to
        % allow r in the input)
        str = str(2:end-1);
        % split the list
        out = strsplit(str,','); % cell array
        num_var = length(out);

        % after splitting string by , split the split strings
        out2 = cellfun(@(x) strsplit(x,'r'),out,'UniformOutput',false); % cell of cells
        nan_var = zeros(num_var,1);
        if_var = cellfun(@length,out2); % length of cell inside cell, 1 or 2
        lb = zeros(nnz(if_var-1),1); % subtract 1 to make values 0 or 1
        ub = zeros(nnz(if_var-1),1);

        a_ = 1;
        for a = 1:num_var
           if (if_var(a)-1) == 1
               % if a variable range is specified (i.e. a -> b)
               nan_var(a) = NaN;
               lb(a_) = str2double(out2{a}{1});
               ub(a_) = str2double(out2{a}{2});
               a_ = a_ + 1;
           else
               % else if the variable is to be fixed
               nan_var(a) = str2double(out2{a});
           end
        end
    else
        % don't process the string
        nan_var = str;
        lb = [];
        ub = [];
    end
    
    % output in rows to be consistent with input
    nan_var = reshape(nan_var,1,[]);
    lb = reshape(lb,1,[]);
    ub = reshape(ub,1,[]);
end