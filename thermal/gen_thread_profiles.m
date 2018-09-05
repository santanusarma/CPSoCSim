% GEN_THREAD_PROFILES generates random thread power and performance
% profiles that can then be used to simulate task migration policies. 
%
% GEN_THREAD_PROFILES requires the following global variables to be defined
% in the calling function:
%   block_power_file - specifies the relative power consumption of the
%       functional blocks in each core (see task_migration_explore.m)
%   P_dyn_max, P_sta_max - mean of the total dynamic and static power per
%       core
%   N_core_x, N_core_y - number of cores in each direction
%
%   thread_info = GEN_THREAD_PROFILES(N_thread, IPC_mean, 
%   interval_mean_ratio_IPC, interval_mean_ratio_P) takes as input:
%       N_thread - number of threads,
%       IPC_mean - mean IPC over all tasks
%       interval_mean_ratio_IPC - the ratio of the interval and the mean
%           for the uniform random distibution used to generate thread IPC
%       interval_mean_ratio_P - - the ratio of the interval and the mean
%           for the uniform random distibution used to generate thread
%           power profiles
%    and returns a structure called thread_info with the following fields:
%           IPC - an array of length N_thread, with the i^th element
%               indicating the IPC of the i^th thread.
%           initialCoreRow - an integer array of length N_thread, with the
%               i^th element indicating the core row to which the i^th
%               thread is initially assigned.
%           initialCoreColumn - an integer array of length N_thread, with the
%               i^th element indicating the core column to which the i^th
%               thread is initially assigned.
%           P_dyn - a matrix of size N_units_per_core x N_thread, with the
%               i^th column listing the dynamic power consumed by each 
%               functional block of any core that runs the i^th thread at 
%               full speed.
%           P_sta - a matrix of size N_units_per_core x N_thread, with the
%               i^th column listing the static power consumed by each 
%               functional block of any core that runs the i^th thread.
%
% Ravishankar Rao, Arizona State University
% Created, Oct 22, 2007

function thread_info = gen_thread_profiles(N_thread, IPC_mean, interval_mean_ratio_IPC, interval_mean_ratio_P)

global block_power_file P_dyn_max P_sta_max N_core_x N_core_y;

[label, p_dyn, p_sta, p_stb] = textread(block_power_file, '%s\t%f\t%f\t%f');
N_units_per_core = length(p_dyn);

r = 1-interval_mean_ratio_IPC/2 + interval_mean_ratio_IPC*rand(N_thread, 1);
IPC = IPC_mean*r;
r = 1-interval_mean_ratio_P/2 + interval_mean_ratio_P*rand(N_units_per_core, N_thread);
P_dyn = P_dyn_max*r.*repmat(p_dyn, 1, N_thread);
r = 1-interval_mean_ratio_P/2 + interval_mean_ratio_P*rand(N_units_per_core, N_thread);
P_sta = P_sta_max*r.*repmat(p_sta, 1, N_thread);

initialCoreRow = zeros(N_thread, 1);
initialCoreColumn = zeros(N_thread, 1);
i = 1;
while i <= N_thread
    % row ~ U[1, N_core_y], col ~ U[1, N_core_x]; row, col are integers
    row = 1 + floor(N_core_y*0.99*rand(1, 1));
    col = 1 + floor(N_core_x*0.99*rand(1, 1));
    initialCoreRow(i) = row;
    initialCoreColumn(i) = col;
    
    % Check if this random combination occured before
    duplicate = 0;
    if i > 1
        for j = 1:(i-1)
            if (row == initialCoreRow(j)) && (col == initialCoreColumn(j))
                duplicate = 1;
                break;
            end
        end
    end
    if duplicate ~= 1
        i = i + 1;
    end
end

thread_info = struct('IPC', IPC, 'initialCoreRow', initialCoreRow, 'initialCoreColumn', initialCoreColumn, 'P_dyn', P_dyn, 'P_sta', P_sta);

    