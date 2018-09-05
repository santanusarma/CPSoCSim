% GET_POWER_TEMP computes the power consumption and the resulting
% temperature of every functional block for the corresponding input speed
% of cores.
%
%   [T,P,vdd] = get_power_temp(spd, Sim_params, Power)
%
%   reads simulation and model parameters from 'Sim_params', and power
%   profile from 'Power' structures respectively. To learn more about the 
%   structures, type'help collect_sim_params.m'. It also reads the input speed
%   'spd'. It also obtains the initial temperature though 'Sim_params'.
%
%   It is capable of computing power for both 'DVFS' and migration. 
%   It also computes the power for only one core or for all cores. The options
%   for these are set in Sim_params.
%
%   The output consists of temperature, power and voltage.
%
% Vinay Hanumaiah, Arizona State University

function [T,P,vdd] = get_power_temp(spd, Sim_params, Power)
    core_units = Sim_params.core_units;
    comp_method = Power.comp_method;
    n_core = Sim_params.n_core;
    T_max = Sim_params.T_max;
    
    if strcmp(Sim_params.dtm_op, 'migration')
        comp_method = 'core';
        curr_task = Sim_params.curr_task;

        T_prev = Sim_params.T_max*ones(Sim_params.core_units,1);
        vdd = interp2(Power.T_range, Power.spd_range, Power.vdd_tab, ...
            max(T_prev), min(spd,1), 'linear', 1);

        P = get_total_power();
        T = Sim_params.G_app_inv *(Sim_params.P_app' + P);
    % get total power including leakage power
        P = P + diag(Sim_params.G_ldt)*T; 
    elseif strcmp(Sim_params.dtm_op, 'scaling')
        T_prev = Sim_params.scaling.T;

        if strcmp(comp_method, 'core')
            tmpr = max(T_prev);
            curr_task = Sim_params.curr_core;
            R = Sim_params.scaling.R(:,:,curr_task);
        elseif strcmp(comp_method, 'total')
            tmpr = zeros(Sim_params.n_core,1);
            for core = 1:Sim_params.n_core
                units_range = (core-1)*core_units+1:core*core_units;
                tmpr(core) = max(T_prev(units_range));
            end
            R = Sim_params.R;
        end
        
        vdd = interp2(Power.T_range, Power.spd_range, Power.vdd_tab, ...
                min(tmpr,T_max), min(spd,1), 'linear', 1);

        P = get_total_power();
        T = R*P;
    end
    
% sub-function to obtain total power
    function [P] = get_total_power()
        % first obtain dynamic power
        if strcmp(comp_method, 'total')
            P = Power.P_d_tot .* (Power.D * (spd.*vdd.^2));
            P_range = (1:n_core*core_units)';
        elseif strcmp(comp_method, 'core') || strcmp(comp_method, 'mig-total')
            P = Power.P_d(:,curr_task)*(spd*vdd^2);
            P_range = (1:core_units)';
        end
    
        % add leakage power obtained through interpolation of
        % 'speed-voltage-temperature' table
        if strcmp(Sim_params.dtm_op, 'scaling')
            P(P_range) = P(P_range) + interp3(Power.T_range, ...
                (1:n_core*core_units)', Power.Vdd_range, Power.P_s, ...
                min(T_prev(P_range),T_max), P_range, ...
                reshape(repmat(vdd,core_units,1),length(P_range),1), ...
                'linear', 0.15);
        elseif strcmp(Sim_params.dtm_op, 'migration')
           P(P_range) = P(P_range) + interp3(Power.T_range, ...
                (1:n_core*core_units)', Power.Vdd_range, ....
                Power.P_s(1:n_core*core_units,:,:), zeros(length(P_range),1), ...
                P_range, reshape(repmat(vdd,core_units,1),length(P_range),1), ...
                'linear', 0.15);
        end
     end
end
