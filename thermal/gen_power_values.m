function [P_d, P_s] = gen_power_values()


%==========================================================================
%ADDITIONAL HELP & EXAMPLE USAGE
%==========================================================================
%Version    : 1.0
%Date       : July 2012
%Author     : Santanu Sarma, University of California, Irvine
%Address    : 3069, Embedded System Lab, Bren Hall, UCI, CA 92697
%Email      : santanus@uci.edu
%Web        : https://students.ics.uci.edu/~santanus/
%==========================================================================


% This functions compute the dynamic power and subthreshold leakage for the
% Alpha EV6 processor. The output is restuned in P_d and P_s. The values are
% on a mat file Power.mat 

% INPUTS
% Nil 
% This should have been customized for multi_core floor plans

% OUTPUTS
% P_d = dynamic power, and 1d array of size of the number of blocks in the
% processor
% P_s = subthreshold power, 3d array.  The size is a function of no of 
% functionsla units, Volatge and temeperature range
% Power.mat save to a mat file 

% EXAMPLE
% Need major modification or update 
% see example usage gen_power_values_test.m
    %% Enum modules
    Constant = struct(...
        'MODULE_DECODE', 1, ...
        'MODULE_BRANCH', 2, ...
        'MODULE_RAT', 3, ...
        'MODULE_RUU', 4, ...
        'MODULE_LSQ', 5, ...
        'MODULE_IALU1', 6, ...
        'MODULE_IALU2', 7, ...
        'MODULE_IALU3', 8, ...
        'MODULE_IALU4', 9, ...
        'MODULE_FALU1', 10, ...
        'MODULE_FALU2', 11, ...
        'MODULE_INTREG', 12, ...
        'MODULE_FPREG', 13, ...
        'MODULE_ITLB', 14, ...
        'MODULE_IL1', 15, ...
        'MODULE_DTLB', 16, ...
        'MODULE_DL1', 17, ...
        'PATH_LENGTH', 256, ...
        'LINE_LENGTH', 256, ...
        'TERM_LENGTH', 64, ...
        'DATA_WIDTH', 64, ...
        'ADDR_WIDTH', 32, ...
        'DCL_GATE_COUNT', 160, ...
        'RUU_SELECTION_GATE_COUNT', 189, ...
        'IALU_GATE_COUNT', 68215, ...
        'FALU_GATE_COUNT', 128369, ...
        'DYNAMIC_SWITCH_FACTOR', 0.15, ...      %
        'C_GATE', (0.75 + 0.91) * 1e-15,  ...   % gate cap as 0.91 fF
        'GATE_AREA', 6.5e-12,  ...              % 6.5 um^2
        'SRAM_CELL_AREA', 2e-12,  ...           % 2 um^2
        'MEM_DENSITY', 0.8, ...
        'DATA_DENSITY', 1.0, ...    
        'MD_NUM_IREGS',	32, ...
        'MD_NUM_FREGS',	32, ...
        'CACHE_POWER_FACTOR', 3.0, ...
        'RUU_size', 8, ...
        'LSQ_size', 8, ...
        'APPEND', 0, ...
        'ASSIGN', 1, ...
        'SF_DYN', 3.5, ...                      % scale factor dynamic power
        'SF_STD', 6.0, ...                      % scale factor subthreshold ??
        'opcode_length', 32, ...
        'gate_count', 68215, ...
        'sram_pd_a', 5.51746e-10, ...
        'sram_pd_b', 9.31735e-8, ...
        'sram_ps_a', -614.9807, ...
        'sram_ps_b', 3528.4329, ...
        'sram_ps_x', 5.2972e-9, ...
        'sram_ps_y', 1.2080e-8, ...
        'sram_ps_z', 5.2946e-9, ...
        'sram_ps_c', -711.9226, ...
        'sram_ps_d', 3725.5342, ...
        'Total_modules', 20);

    K = 273.15;
    
    Params = struct(...
        'opt', 'dynamic', ...
        'tmpr', 300, ...
        'vdd', 1);
    
    Rent_params = struct(...
    'LOCAL_INTERMEDIATE_BOUNDARY', 21, ...
    'CW_LOCAL', 176.188e-12, ...
    'CW_INT', 180.068e-12, ...
    'ALPHA', 1.0, ...
    'K', 4.0, ...
    'P', 0.55);
    
    %% Dynamic power values
    Params.type = 'dynamic';
    Params.Rent_params = Rent_params;
    P_d = module_list_comp(Constant, Params);
    
    %% Leakage power values
    Params.type = 'leakage';
    vdd_range = 0.7:0.05:1.4;
    tmpr_range = 30+K:5:110+K;
    vdd_length = length(vdd_range);
    tmpr_length = length(tmpr_range);
    
    P_s = zeros(vdd_length, tmpr_length, Constant.Total_modules);

    for vdd_index = vdd_length
        for tmpr_index = tmpr_length
            Params.vdd = vdd_range(vdd_index);
            Params.tmpr = tmpr_range(tmpr_index);
            P_s(vdd_index, tmpr_index, :) = module_list_comp(Constant, Params);
         end
    end

    save Power P_d P_s;
    
    
    
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

    function [Power] = module_list_comp(Constant, Params)
        Power = zeros(Constant.Total_modules, 1);
        %IALUs
        Params.gate_count = Constant.IALU_GATE_COUNT;
        for module = Constant.MODULE_IALU1:Constant.MODULE_IALU4 
            Power(module) = non_memory_power(Constant, Params);
        end
        
        % FALUs
        Params.gate_count = Constant.FALU_GATE_COUNT;
        for module = Constant.MODULE_FALU1:Constant.MODULE_FALU2
            Power(module) = non_memory_power(Constant, Params);
        end
        % rename table
        Power(Constant.MODULE_RAT) = rename_table_power(Constant, Params);
        
        % DECODE
        Params.prev_power = 0;
        Params.row = 65;
        Params.col = Constant.opcode_length;
        Params.opt = Constant.ASSIGN;
        % shouldn't be the issue_width 
        Power(Constant.MODULE_DECODE) = sram_array_power(Constant, Params); 
        
        % INT_REG
        Params.prev_power = 0;
        Params.row = Constant.MD_NUM_IREGS;
        Params.col = Constant.DATA_WIDTH;
        Params.opt = Constant.ASSIGN;
        Power(Constant.MODULE_INTREG) = sram_array_power(Constant, Params);
        
        % FP_REG
        % Other parameters remain same as INT_REG
        Params.row = Constant.MD_NUM_FREGS; 
        Power(Constant.MODULE_FPREG) = sram_array_power(Constant, Params);
        
        % RUU_SELECTION
        Params.gate_count = Constant.RUU_SELECTION_GATE_COUNT;
        Params.prev_power = non_memory_power(Constant, Params);        
        % RUU_WAKEUP
        Params.row = Constant.RUU_size;
        Params.col = ceil(log2(Params.row));
        Params.opt = Constant.APPEND;
        sram_array_power(Constant, Params);
        % RUU_RS
        Params.col = Constant.DATA_WIDTH;
        Power(Constant.MODULE_RUU) = sram_array_power(Constant, Params);
        
        % LSQ_WAKEUP. CAM is also modeled by sram array
        % FIXME: address in LSQ!
        Params.prev_power = 0;
        Params.row = Constant.LSQ_size;
        Params.col = Constant.DATA_WIDTH;
        Params.opt = Constant.ASSIGN;
        Params.prev_power = sram_array_power(Constant, Params);
        % LSQ_RS 
        Params.opt = Constant.APPEND;
        Power(Constant.MODULE_LSQ) = sram_array_power(Constant, Params);
    end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

    function [Power] = rename_table_power(Constant, Params)
        % Dependence Check Logic
        Params.gate_count = Constant.DCL_GATE_COUNT;
        Power = non_memory_power(Constant, Params);
        % RAT
        Params.prev_power = Power;
        Params.row = Constant.MD_NUM_IREGS;
        Params.col = ceil(log2(Params.row));
        Params.opt = Constant.APPEND;
        sram_array_power(Constant, Params);
    end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

    function [Power] = non_memory_power(Constant, Params)
        SF_DYN = Constant.SF_DYN;
        SF_STD = Constant.SF_STD;
        sf = Constant.DYNAMIC_SWITCH_FACTOR;
        c_gate = Constant.C_GATE;
        gate_count = Params.gate_count;
        tmpr = Params.tmpr;
        vdd = Params.vdd;
        
        switch Params.type
            case 'dynamic'
                gate_power = SF_DYN * (0.5 * gate_count * sf * c_gate);
                interconnect_power = SF_DYN * 0.5 * sf * ...
                    rents(Constant, Params.Rent_params, gate_count);

                Power = (gate_power + interconnect_power);

            case 'leakage'
                Power = SF_STD * gate_count * vdd * (2.48325558 + ...
                    2.26321378 + 2.6529006) / 3.0 * 1e-8 * tmpr^2 * ...
                    exp(-(Constant.sram_ps_a * vdd + ...
                    Constant.sram_ps_b) / tmpr);
        end
    end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

    function [Power] = sram_array_power(Constant, Params)
        opt = Params.opt;
        prev_power = Params.prev_power;
        
    	% default is for ASSIGN
        append_factor = 0.0;
        scale_factor = 1.0;
    
        if opt == Constant.ASSIGN
            append_factor = 0.0;
            scale_factor = 1.0; 
        elseif opt == Constant.APPEND
            append_factor = 1.0;
            scale_factor = 1.0; 
        elseif opt > 1
            error('opt > 1');
        else
            append_factor = 0.0;
            scale_factor = opt;
        end
          
        switch Params.type
            case 'dynamic'
                Power = prev_power * append_factor + (1 + scale_factor) * ...
                        power_array(Constant, Params);
            case 'leakage'
                Power = prev_power * append_factor + (1 + scale_factor) * ...
                        power_array(Constant, Params);
        end
    end
%-------------------------------------------------------------------------- 
%--------------------------------------------------------------------------

    function [Power] = power_array(Constant, Params)
        SF_STD = Constant.SF_STD;
        tmpr = Params.tmpr;
        vdd = Params.vdd;
        row = Params.row;
        col = Params.col;

        switch Params.type
            case 'dynamic'
                  Power = Constant.SF_DYN * (Constant.sram_pd_a * row + ...
                    Constant.sram_pd_b * col) * Constant.CACHE_POWER_FACTOR;
            case 'leakage'
                std_power_circuits = SF_STD * vdd * (Constant.sram_ps_x * ...
                    row * col + Constant.sram_ps_y * col) * tmpr^2 * ...
                    exp(-(Constant.sram_ps_a * vdd + Constant.sram_ps_b) / tmpr);
                std_power_cells = SF_STD * vdd * (Constant.sram_ps_z * ...
                    row * col) * tmpr^2 * exp(-(Constant.sram_ps_c * vdd + ...
                    Constant.sram_ps_d) / tmpr);
    
                Power = (std_power_circuits + std_power_cells) * ...
                    Constant.CACHE_POWER_FACTOR; 
        end
    end
%-------------------------------------------------------------------------- 
%--------------------------------------------------------------------------

    function [total_cap] = rents(Constant, Rent_params, gate_count)
        gate_pitch = 2.0 * sqrt(Constant.GATE_AREA);
        total_cap = 0.0;
    
        for n = 1.0:2.0 * sqrt(gate_count)
            if n >= Rent_params.LOCAL_INTERMEDIATE_BOUNDARY
                break;
            end
            idf = get_idf(Rent_params, n, gate_count);
            total_cap = total_cap + idf * (n * gate_pitch * ...
                        Rent_params.CW_LOCAL);
        end
        
        for n = n:2 * sqrt(gate_count)
            idf = get_idf(Rent_params, n, gate_count);
            total_cap = total_cap + idf * (n * gate_pitch * ...
                        Rent_params.CW_INT);
        end
    end
 
%--------------------------------------------------------------------------
% This function computes idf 
% INPUTS:
% Rent_params : a structue conating  fields, purpose ??
% L = ?
% N = ?
% OUTPUTS:
% Idf =? some sort of leakage current , units mA?
% EXAMPLE:
%--------------------------------------------------------------------------

    function [idf] = get_idf(Rent_params, L, N)
        ALPHA = Rent_params.ALPHA;
        P = Rent_params.P;
        K = Rent_params.K;
        
        tau = get_tau(Rent_params, N);
    
        if L < sqrt(N)
            a = ALPHA * K / 2;
            b = tau * (L^3 /3 - 2*sqrt(N) * L^2 + 2*N*L);
            c = L^(2*P-4);
            idf = a*b*c;
        else
            a = ALPHA * K / 6;
            b = tau * (2*sqrt(N) - L)^3;
            c = L^(2*P-4);
            idf = a*b*c;
        end
    end
    

%--------------------------------------------------------------------------
%This function computes tau. The expression of tau is used in the current
%computation 
%
% INPUTS :
%  Rent_params      = parameters for ??
%  N                = ??

% OUTOUTS:
% tau               = a constant used in idf computation ??
%--------------------------------------------------------------------------

    function [tau] = get_tau(Rent_params, N)
        
        P = Rent_params.P;
        
        if P == 0.5
            a = 4 * N - 4 * sqrt(N);
            b = sqrt(N) * (-2 * log(N) - 6 + 2 * log(4)) + 4*N - 2/3;
            tau = a/b;
        else
            a = 2 * N * (1 - N^(P-1));
            b = - N^P * (1 + 2*P - 2^(2*P-1)) / (P*(2*P-1) * (P-2) * (2*P-3));
            c = - 1 / (6*P) + 2 * sqrt(N) / 2^(2*P-1) - N / (P-1);
            tau = a / (b + c);
        end
    end
    
end
