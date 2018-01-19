function [param,lb,ub] = mcVary(param,vars)
% uses split_input, another function, to separate a r b into range of a - b

    if isempty(vars)
        lb = cell(1,21);
        ub = cell(1,21);
        % splits the 'r' strings into NaN and lb/ub ranges
        [param.comp_x,lb{1},ub{1}] = split_input(param.comp_x);
        [param.comp_y,lb{2},ub{2}] = split_input(param.comp_y);
        [param.comp_z,lb{3},ub{3}] = split_input(param.comp_z);
        [param.imp_x,lb{4},ub{4}] = split_input(param.imp_z);
        [param.imp_y,lb{5},ub{5}] = split_input(param.imp_y);
        [param.imp_z,lb{6},ub{6}] = split_input(param.imp_z);
        [param.Ex,lb{7},ub{7}] = split_input(param.Ex);
        [param.Ey,lb{8},ub{8}] = split_input(param.Ey);
        [param.Ez,lb{9},ub{9}] = split_input(param.Ez);
        [param.kick,lb{10},ub{10}] = split_input(param.kick);
        [param.photon,lb{11},ub{11}] = split_input(param.photon);
        [param.comp_range,lb{12},ub{12}] = split_input(param.comp_range);
        [param.imp_range,lb{13},ub{13}] = split_input(param.imp_range);
        [param.Lx,lb{14},ub{14}] = split_input(param.Lx);
        [param.Ly,lb{15},ub{15}] = split_input(param.Ly);
        [param.Lz,lb{16},ub{16}] = split_input(param.Lz);
        [param.T0,lb{17},ub{17}] = split_input(param.T0);
        [param.Te,lb{18},ub{18}] = split_input(param.Te);
        [param.ni,lb{19},ub{19}] = split_input(param.ni);
        [param.co,lb{20},ub{20}] = split_input(param.co);
        [param.specularity,lb{21},ub{21}] = split_input(param.specularity);
        lb = cell2mat(lb);
        ub = cell2mat(ub);
    else
        lb = [];
        ub = [];
        % assign the variables to the NaN
        var_counter = cumsum([sum(isnan(param.comp_x)) sum(isnan(param.comp_y)) sum(isnan(param.comp_z))...
                              sum(isnan(param.imp_x)) sum(isnan(param.imp_y)) sum(isnan(param.imp_z))...
                              sum(isnan(param.Ex)) sum(isnan(param.Ey)) sum(isnan(param.Ez))...
                              sum(isnan(param.kick)) sum(isnan(param.photon)) sum(isnan(param.comp_range)) sum(isnan(param.imp_range))...
                              sum(isnan(param.Lx)) sum(isnan(param.Ly)) sum(isnan(param.Lz)) sum(isnan(param.T0)) sum(isnan(param.Te))...
                              sum(isnan(param.ni)) sum(isnan(param.co)) sum(isnan(param.specularity))]);
        param.comp_x(isnan(param.comp_x)) = vars(1:var_counter(1));
        param.comp_y(isnan(param.comp_y)) = vars(var_counter(1)+1:var_counter(2));
        param.comp_z(isnan(param.comp_z)) = vars(var_counter(2)+1:var_counter(3));
        param.imp_x(isnan(param.imp_x)) = vars(var_counter(3)+1:var_counter(4));
        param.imp_y(isnan(param.imp_y)) = vars(var_counter(4)+1:var_counter(5));
        param.imp_z(isnan(param.imp_z)) = vars(var_counter(5)+1:var_counter(6));
        param.Ex(isnan(param.Ex)) = vars(var_counter(6)+1:var_counter(7));
        param.Ey(isnan(param.Ey)) = vars(var_counter(7)+1:var_counter(8));
        param.Ez(isnan(param.Ez)) = vars(var_counter(8)+1:var_counter(9));
        param.kick(isnan(param.kick)) = vars(var_counter(9)+1:var_counter(10));
        param.photon(isnan(param.photon)) = vars(var_counter(10)+1:var_counter(11));
        param.comp_range(isnan(param.comp_range)) = vars(var_counter(11)+1:var_counter(12));
        param.imp_range(isnan(param.imp_range)) = vars(var_counter(12)+1:var_counter(13));
        param.Lx(isnan(param.Lx)) = vars(var_counter(13)+1:var_counter(14));
        param.Ly(isnan(param.Ly)) = vars(var_counter(14)+1:var_counter(15));
        param.Lz(isnan(param.Lz)) = vars(var_counter(15)+1:var_counter(16));
        param.T0(isnan(param.T0)) = vars(var_counter(16)+1:var_counter(17));
        param.Te(isnan(param.Te)) = vars(var_counter(17)+1:var_counter(18));
        param.ni(isnan(param.ni)) = vars(var_counter(18)+1:var_counter(19));
        param.co(isnan(param.co)) = vars(var_counter(19)+1:var_counter(20));
        param.specularity(isnan(param.specularity)) = vars(var_counter(20)+1:var_counter(21));
    end
end