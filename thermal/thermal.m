%THERMAL SIMULATION of Alpha Processor (1 core only)
% EXTRACT_HOTSPOT_MODEL runs Hotspot on a given floorplan file, and 
% extracts the resulting thermal RC circuit in state space form:
%   dT/dt = A_HS T + B_HS P
% where T and P are vectors of length = (number of blocks in floorplan +
% number of blocks in Hotspot package).
%% Extract the A and B matrices for the state space representation of the
% RC thermal equivalent circuit. Hotspot uses the following state-space
% equation: dT/dt + inv(A) * B * T = inv(A) * P, 
% or dT/dt = -inv(A) * B * T + inv(A) * P. Comparing with the standard
% state space notation: dT/dt = A_HS * T + B_HS * P, we have
% A_HS = -inv(A) * B, B_HS = inv(A). The inv(A) and B matrices are computed
% by Hotspot and dumped into the files matrix_inva and matrix_b,
% respectively. 

% 
% %Load the thermal model and the power values
% load('/Users/santanusarma/Dropbox/virtual sensor/MATLAB/magma_v2/models/1x1/ev6_1x1.mat')
% load('Power.mat')
% 
% Thermal_model
% A_HS=Thermal_model.A_HS
% B_HS=Thermal_model.B_HS
% [ax ay]=size(A_HS)
% [bx by]=size(B_HS)
% %define the output equation, measurement matric C_HS, D_HS
% % The measurement matrix 
% C_HS=eye(ax);%ones(1,ax)
% % The input matrix should only contain element equal to no of power gen
% % blocks
% D_HS=eye(ax)
% 
% % The size of the matrices are A=54x54, B is 54x54
% % This size is achievd as for 20 FU's ie 2*FU+14= 20*2+ 14
% 
% 
% 
% %Create a statespace model
% SysAlpha=ss(A_HS, B_HS, C_HS, D_HS);
% 
% %Convert the continuous time model to discrete model 
% Ts=1 ; % Sampling time
% SysAlphaD=c2d(SysAlpha, Ts);
% %Check The step response
% SimTime=100; %Sec
% %hs1=tf(A_HS, B_HS, C_HS, D_HS, 1)
% %step(SysAlpha)
% 
% 
% 
% %Visualize the power for FUs
% %==========================================================================
% 
% %Check all in the included function in the simulator and update the code
% 
% % Test 1: (Pass)
% % leakage_surface
% % This script plots the leakage power surface
% leakage_surface
% 
% % Test 2:
% 
% 
% 
% %Visualize the Floorplan
% %==========================================================================
% % Test 1: (Pass)
% 
% % Test for create_spatial_plot
% % Input is temperature , Output is a surface plot 
% T=30;
% plot_mat = create_spatial_plot(T)
% 
% 
% % Test 2: (Pass)
% % create_floorplan
%  create_floorplan ('ev6_1x1.flp', 2,2, 'mc_ev6_2x2.flp')
%  
%  
%  %Test 3: Floor Plan Viewing (Pass)
%  % function view_floorplan(w, h, x, y, hFigure, label, z)
%  view_floorplan_test;
%  
%  
%  pause(2)
%  close all;
 
%==========================================================================
% Simulation parameters
%==========================================================================
clear, clc

% Params.n_core_x =2;
% Params.n_core_y =2
% Params.thread_map=[1 2 3 4]
% Params.n_task =4
% Params.f_max = 2e9
% Params.vdd_max = 1.2
% Params.T_ini =0
% Params.IPC = [1.8; 1.85; 1.95; 2.45]
% Params.t_step=2
% Params.policy = 'zero_slack'
% Params.dtm_op = 'scaling'
% Params.method = 'dvfs'
% Params.sub_method = 'accurate'
% Params.exp = 'temporal'
% Params.inst= [0.5; 0.3; 0.7; 0.9]*1e12
% Params.n_states =10
% Params.t_sim =200;
% Params.discrete =false


Params = struct(...
    'n_core_x', 2, ...
    'n_core_y', 2, ...
    'thread_map', [1 2 3 4], ...
    'n_task', 4, ...
    'f_max', 2e9, ...
    'vdd_max', 1.2, ...
    'T_ini', 0, ...
    'IPC', [1.8; 1.85; 1.95; 2.45], ...
    't_step', 2, ...
    'policy', 'zero_slack', ...
    'dtm_op', 'scaling', ...
    'method', 'dvfs', ...
    'sub_method', 'accurate', ...
    'exp', 'temporal', ...
    'inst', [0.5; 0.3; 0.7; 0.9]*1e12, ...
    'n_states', 10,...
    't_sim', 200, ...
    'discrete', false);

n_core = Params.n_core_x*Params.n_core_y;
Params.P_d_sum = 230;%(1+0.2)^(n_core/2)*200;
Params.P_a_sum = 60;


%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================

%Constant = read_model_parameters

%% Die parameters

% Number of cores in multicore processor in horizontal (X) direction
n_core_x = 1; 
% Number of cores in multicore processor in vertical (Y) direction
n_core_y = 1; 
% Should lateral thermal resistances be omitted in thermal model
omit_lateral = false; 
% Should heat spreader and heat sink center be modeled as a single block
flat_pkg_center = true;
% Floorplan name: built-in examples: simple, core_cache, ev6
floorplan_name = 'ev6';
% Use floorplan plan that was saved before instead of creating a new one
use_saved_floorplan = true;

% Minimum (normalized) speed per core, normalized units i.e no unit
s_min = 0.0; % Maximum is 1. Recommended value: 0.1

%% Thermal parameters

% Internal ambient temperature in degree C
t_amb_abs = 35; 
% Maximum allowed die temperature in degree C
t_max_abs = 110; 
% Maximum allowed die temperature in C relative to internal ambient
t_max = t_max_abs - t_amb_abs;
% 0 C in K
t_0_C_in_K = 273;
% Hotspot convection thermal resistance, 0.35 K/W default
r_convec = 0.35;

% Mobility degradation factor (MDF): 
%   Speed <= (Speed at t_amb_abs) * ((t_amb_abs in K) / ...
%           (max. chip. temp. in K))^MDF
% Set MDF to 0 to ignore mobility degradation with temperature

MDF = 0.0; 

% Use previously saved thermal model instead of generating new one 
use_saved_thermal_model = true;

% Should thermal model combine convection and heat sink (default Hotspot 4
% behavior) or should they be modeled as separate nodes
combine_hs_conv = true;

%--------------------------------------------------------------------------
%% Power parameters
%--------------------------------------------------------------------------
% Power data file. Contains reference power consumption values (dynamic, 
% active leakage and standby leakage) of each functional unit for a single 
% core. The dynamic powers are measured at maximum speed, while the leakage
% powers are measured at maximum temperature. These values are then 
% normalized and scaled respectively by p_d_at_fmax_for_SC, 
% p_a_at_Tmax_for_SC, and p_s_at_Tmax_for_SC. To see the file format for 
% the power_data_file, type 'help read_power_data'.

power_data_file = 'power.txt'; % The file containing power values for FU

% Format of the power.txt file 


% max dynamic power in W at max freqeuncy for scaling
p_d_at_fmax_for_SC = 230;  
% max active leakage power in W at max temperature in degree for scaling
p_a_at_Tmax_for_SC = 160;
% max standby leakage power in W at max temperature in degree for scaling
p_s_at_Tmax_for_SC = 10;

% The nonlinear relation between leakage and temperature can be modeled as
% a 'constant' (essentially, no dependence), 'linear', or 'exponential'. 

LDT_model = 'exponential'; % lekage dependent temperature model 


n_PWL_LDT = 5;  %            
n_iter_LP_PWL_LDT = 10; %

use_SPEC_benchmarks = false;
use_MiBench_benchmarks = true;

power_dir = 'power_data';

%--------------------------------------------------------------------------
%Thread Performance Model
%--------------------------------------------------------------------------
%% Thread parameters

% Number of threads. Must be an integer from 1 to n_core_x * n_core_y - 1.
n_thread = n_core_x * n_core_y - 1;
% The Instructions Per Clock (IPC) of each thread can vary spatially (i.e.
% between threads) and temporally (for each thread). These variations are
% modeled as samples drawn from a uniform random distribution. The
% parameters of this distribution are:
%   IPC_mean, and IMR_IPC_spatial (IMR_IPC_temporal) for the spatial
%   variations
%   IMR refers to Interval-Mean Ratio. It is the ratio between the
%   interval of the uniform random distribution and the mean.
% The samples of the thread IPCs are drawn from a uniform random
% distribution over the interval IPC_mean * [1 - IMR_IPC_spatial/2, 1 +
% IMPR_IPC_spatial/2]. Each of these samples then represent the mean taken
% over time for the respective thread. Let this mean IPC for the i^th
% thread be IPC(i). Then, the actual IPC for thread i in each time slot is
% chosen from a uniform distribution over the interval IPC(i) * [1 -
% IMR_IPC_temporal/2, 1 + IMR_IPC_temporal/2]. 
% Note: Currently, IPC variation is not implemented. So, keep
% IMR_IPC_spatial and IMR_IPC_temporal at 0.

IPC_mean = 1.0;         % normalized IPC
IMR_IPC_spatial = 0;    % The random model is not considered
IMR_IPC_temporal = 0;   % The random model is not considered

% SPEC benchmarks
% List of benchmarks whose power profiles are available 
%   1      bzip2
%   2      crafty
%   3      galgel
%   4      gap
%   5      gcc
%   6      mcf
%   7      twolf
% This list corresponds to the files present in the power_data directory.
% If this changes, update the above list, which is sorted in alphabetical
% order. To select bzip2, crafty, mcf, and twolf for 4 cores, set the
% vector Thread_initial_map to [1 2 6 7]'. Threads can be repeated on
% different cores. So, for 9 cores, you could do [1:7 2 3]'.
Thread_initial_map = [1 2 5 7]';
% Ignore this for now
if use_SPEC_benchmarks
    IPC_SPEC = [20892 16210 15230 20211 11139 25630 30453 20892 20892]'/20000;
elseif use_MiBench_benchmarks
    IPC_SPEC = [1.8 1.85 1.95 2.45]';
end

n_benchmarks = length(IPC_SPEC);


% The power values are also varied around the mean values obtained above
% using a random distribution. The parameters and perturbation methods are
% similar to the one mentioned above for IPC.
IMR_power_spatial = 0.0;
IMR_power_temporal = 0.0;

% The random generation algorithm used. Type 'help rand' to see options.
rand_gen_method = 'state';

% The state (or seed) for the random number generator. This ensures
% repeatable experiments.
state_IPC_spatial = 100;
state_IPC_temporal = 200;
state_power_spatial = 301;
state_power_temporal = 400;
state_thread_location = 500;

%% Simulation parameters

% The time interval at which thread migration is performed for time-based
% policies like rotation (round-robin). Must be an integer multiple of
% n_core * mig_overhead.
mig_interval = 10000e-6; % sec
% The time overhead for migrating each thread
mig_overhead = 10e-6; % sec 
% The number of time steps, each of length time_step, to simulate 
n_step = 15000; % Units of each step ??
% Number of initial steps to ignore in computing throughput to remove
% transient effects. Set to zero if none should be ignored.
ignore_first_n_steps = 1;
% Global speed control time window
n_step_per_global_opt = 5000;

% Should the temperature, power, speed, core activity in each time slot be
% saved in the Results data structure
save_profile = true; % Save to memory
% Save to file profile.dat. Use this option when n_step > 10000
save_profile_to_file = true;
% For faster simulation of thread migration policies
ignore_mig_overhead = false; 

%--------------------------------------------------------------------------
% Solver : Lenear Programing (only applicable when MDF =0)
%--------------------------------------------------------------------------
% Choose between accurate lin_prog (linear programming) and almost as
% accurate bin_search (binary search) methods to perform distributed speed
% control for multiple cores. Note that the linear program method is
% currently insensitive to mobility degradation with temperature. So, only
% use lin_prog if MDF = 0.

speed_opt_method = 'analytical';
n_iter_LP_Exp_LDT = 10;
t_tol_LP_Exp_LDT = 0.001;

speed_opt_window = 'local';

% Choose whether to constrain the temperature of some or all of the blocks 
% on the die. Since the thermal upper bound is the same for the entire die, 
% if we could identify the hottest functional unit on each core, it is 
% sufficient to constraint only that unit. This will speed up solving for 
% the optimal speed combination. If constrain_k_hottest_units = 0, all die
% funcitonal units are constrained, else, the n_hottest_units (as
% determined from the temperature vector at the beginning of the time step)
% are constrained. The idea is that only a small set of functional units
% tend to be the hotspots in a core, and that this set does not change much
% between time steps. 
constrain_n_hottest_units = 1;
n_hottest_units = 1;

% Parameters for binary search speed optimization method (See the file
% find_speed_binary_search.m for details).
% Tolerance for normalized speed in local search phase
s_tol = 1e-3; 
% Tolerance for speed scaling factor in global search
f_tol = 1e-2; 

%% Policy parameters

% Thread migration policy
thread_migration_policy = 'none';


% %--------------------------------------------------------------------------
% % Simulation Parameters Structure
% %--------------------------------------------------------------------------
% %% Save parameters in structure
Constant = struct(...
    'n_core_x', n_core_x, ...
    'n_core_y', n_core_y, ...
    'omit_lateral', omit_lateral, ...
    'flat_pkg_center', flat_pkg_center, ...
    'floorplan_name', floorplan_name, ...
    't_amb_abs', t_amb_abs, ...
    't_max_abs', t_max_abs, ...
    't_max', t_max, ...
    't_0_C_in_K', t_0_C_in_K, ...
    'r_convec', r_convec, ...
    'combine_hs_conv', combine_hs_conv, ...
    'power_data_file', power_data_file, ...
    'p_d_at_fmax_for_SC', p_d_at_fmax_for_SC, ...
    'p_a_at_Tmax_for_SC', p_a_at_Tmax_for_SC, ...
    'p_s_at_Tmax_for_SC', p_s_at_Tmax_for_SC, ...
    'LDT_model', LDT_model, ...
    'n_thread', n_thread, ...
    'IPC_mean', IPC_mean, ...
    'IMR_IPC_spatial', IMR_IPC_spatial, ...
    'IMR_IPC_temporal', IMR_IPC_temporal, ...
    'IMR_power_spatial', IMR_power_spatial, ...
    'IMR_power_temporal', IMR_power_temporal, ...
    'rand_gen_method', rand_gen_method, ...
    'state_IPC_spatial', state_IPC_spatial, ...
    'state_IPC_temporal', state_IPC_temporal, ...
    'state_power_spatial', state_power_spatial, ...
    'state_power_temporal', state_power_temporal, ...
    'state_thread_location', state_thread_location, ...
    'mig_interval', mig_interval, ...
    'mig_overhead', mig_overhead, ...
    'ignore_mig_overhead', ignore_mig_overhead, ...
    'n_step', n_step, ...
    'ignore_first_n_steps', ignore_first_n_steps, ...
    'save_profile', save_profile, ...
    'save_profile_to_file', save_profile_to_file, ...
    'speed_opt_method', speed_opt_method, ...
    'MDF', MDF, ...
    'use_saved_floorplan', use_saved_floorplan, ...
    'use_saved_thermal_model', use_saved_thermal_model, ...
    'constrain_n_hottest_units', constrain_n_hottest_units, ...
    'n_hottest_units', n_hottest_units, ...
    's_tol', s_tol, ...
    'f_tol', f_tol, ...
    'thread_migration_policy', thread_migration_policy, ...
    's_min', s_min, ...
    'n_step_per_global_opt', n_step_per_global_opt, ...
    'speed_opt_window', speed_opt_window, ...
    'use_SPEC_benchmarks', use_SPEC_benchmarks, ...
    'use_MiBench_benchmarks', use_MiBench_benchmarks, ...
    'IPC_SPEC', IPC_SPEC, ...
    'n_benchmarks', n_benchmarks, ...
    'Thread_initial_map', Thread_initial_map, ...
    'power_dir', power_dir, ...
    'n_iter_LP_Exp_LDT', n_iter_LP_Exp_LDT, ...
	't_tol_LP_Exp_LDT', t_tol_LP_Exp_LDT, ...
	'n_PWL_LDT', n_PWL_LDT, ...
	'n_iter_LP_PWL_LDT', n_iter_LP_PWL_LDT);


% Constant.n_core_x= n_core_x, ...
% Constant.n_core_y= n_core_y, ...
% Constant.omit_lateral= omit_lateral, ...
% Constant.flat_pkg_center=flat_pkg_center, ...
% Constant.floorplan_name= floorplan_name, ...
% Constant.t_amb_abs=t_amb_abs, ...
% Constant.t_max_abs= t_max_abs, ...
% Constant.t_max= t_max, ...
% Constant.t_0_C_in_K= t_0_C_in_K, ...
% Constant.r_convec= r_convec, ...
% Constant.combine_hs_conv= combine_hs_conv, ...
% Constant.power_data_file= power_data_file, ...
% Constant.p_d_at_fmax_for_SC= p_d_at_fmax_for_SC, ...
% Constant.p_a_at_Tmax_for_SC= p_a_at_Tmax_for_SC, ...
% Constant.p_s_at_Tmax_for_SC= p_s_at_Tmax_for_SC, ...
% Constant.LDT_model= LDT_model, ...
% Constant.n_thread= n_thread, ...
% Constant.IPC_mean = IPC_mean, ...
% Constant.IMR_IPC_spatial= IMR_IPC_spatial, ...
% Constant.IMR_IPC_temporal= IMR_IPC_temporal, ...
% Constant.IMR_power_spatial= IMR_power_spatial, ...
% Constant.IMR_power_temporal= IMR_power_temporal, ...
% Constant.rand_gen_method=rand_gen_method, ...
% Constant.state_IPC_spatial= state_IPC_spatial, ...
% Constant.state_IPC_temporal= state_IPC_temporal, ...
% Constant.state_power_spatial= state_power_spatial, ...
% Constant.state_power_temporal=state_power_temporal, ...
% Constant.state_thread_location= state_thread_location, ...
% Constant.mig_interval= mig_interval, ...
% Constant.mig_overhead= mig_overhead, ...
% Constant.ignore_mig_overhead= ignore_mig_overhead, ...
% Constant.n_step= n_step, ...
% Constant.ignore_first_n_steps= ignore_first_n_steps, ...
% Constant.save_profile= save_profile, ...
% Constant.save_profile_to_file=save_profile_to_file, ...
% Constant.speed_opt_method= speed_opt_method, ...
% Constant.MDF= MDF, ...
% Constant.use_saved_floorplan=use_saved_floorplan, ...
% Constant.use_saved_thermal_model=use_saved_thermal_model, ...
% Constant.constrain_n_hottest_units= constrain_n_hottest_units, ...
% Constant.n_hottest_units=n_hottest_units, ...
% Constant.s_tol=s_tol, ...
% Constant.f_tol= f_tol, ...
% Constant.thread_migration_policy=thread_migration_policy, ...
% Constant.s_min= s_min, ...
% Constant.n_step_per_global_opt=n_step_per_global_opt, ...
% Constant.speed_opt_window=speed_opt_window, ...
% Constant.use_SPEC_benchmarks= use_SPEC_benchmarks, ...
% Constant.use_MiBench_benchmarks=use_MiBench_benchmarks, ...
% Constant.IPC_SPEC= IPC_SPEC, ...
% Constant.n_benchmarks= n_benchmarks, ...
% Constant.Thread_initial_map=Thread_initial_map, ...
% Constant.power_dir= power_dir, ...
% Constant.n_iter_LP_Exp_LDT= n_iter_LP_Exp_LDT, ...
% Constant.t_tol_LP_Exp_LDT= t_tol_LP_Exp_LDT, ...
% Constant.n_PWL_LDT=n_PWL_LDT, ...
% Constant.n_iter_LP_PWL_LDT=n_iter_LP_PWL_LDT
% 


%==========================================================================


%--------------------------------------------------------------------------
% Get parameters
vdd_max = Params.vdd_max;
n_task = Params.n_task;
t_step = Params.t_step;
dtm_op = Params.dtm_op;
sub_method = Params.sub_method;

Constant.LDT_model = 'linear';
Constant.thread_migration_policy = 'none';
Constant.ignore_first_n_steps = 1;
Constant.save_profile = 1;
Constant.save_profile_to_file = 0;
Constant.ignore_mig_overhead = true;
Constant.omit_lateral = 0;
Constant.p_d_at_fmax_for_SC = Params.P_d_sum;
Constant.p_a_at_Tmax_for_SC = Params.P_a_sum;
Constant.use_saved_thermal_model = true;
Constant.mig_interval = t_step;

switch dtm_op
    case {'scaling', 'dse'}
        Constant.n_core_x = Params.n_core_x;
        Constant.n_core_y = Params.n_core_y;
        Constant.n_thread = Constant.n_core_x * Constant.n_core_y;
        Constant.Thread_initial_map = Params.thread_map;
        Constant.speed_opt_window = 'local_global';
        Constant.speed_opt_method = 'LP_NLP';
        Constant.max_inst = ceil(max(Params.inst./Constant.IPC_SPEC)/10e6);
        Constant.use_SPEC_benchmarks = true;
        Constant.use_MiBench_benchmarks = true;
    case 'migration'
        Constant.n_core_x = Params.n_core_x;
        Constant.n_core_y = Params.n_core_y;
        Constant.Thread_initial_map = Params.thread_map;
        Constant.speed_opt_window = 'analytical';
        Constant.speed_opt_method = 'local';
        Constant.max_inst = ceil(max(Params.t_sim./Constant.IPC_SPEC)*1e3);
        Constant.use_SPEC_benchmarks = true;
        Constant.use_MiBench_benchmarks = true;
end


%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%--------------------------------------------------------------------------
% [Die, Power1, Thermal, Thread, Simulation, Policy]= ...
%    compute_model_parameters(Constant)
%==========================================================================

%% Read constants specified in read_model_parameters.m
n_core_x = Constant.n_core_x;
n_core_y = Constant.n_core_y;
omit_lateral = Constant.omit_lateral;
flat_pkg_center = Constant.flat_pkg_center;
floorplan_name = Constant.floorplan_name;
t_amb_abs = Constant.t_amb_abs;
t_max_abs = Constant.t_max_abs;
t_max = Constant.t_max;
t_0_C_in_K = Constant.t_0_C_in_K;
r_convec = Constant.r_convec;
combine_hs_conv = Constant.combine_hs_conv;
power_data_file = Constant.power_data_file;
p_d_at_fmax_for_SC = Constant.p_d_at_fmax_for_SC;
p_a_at_Tmax_for_SC = Constant.p_a_at_Tmax_for_SC;
p_s_at_Tmax_for_SC = Constant.p_s_at_Tmax_for_SC;
LDT_model = Constant.LDT_model;
n_thread = Constant.n_thread;
IPC_mean = Constant.IPC_mean;
IMR_IPC_spatial = Constant.IMR_IPC_spatial;
IMR_IPC_temporal = Constant.IMR_IPC_temporal;
IMR_power_spatial = Constant.IMR_power_spatial;
IMR_power_temporal = Constant.IMR_power_temporal;
rand_gen_method = Constant.rand_gen_method;
state_IPC_spatial = Constant.state_IPC_spatial;
state_IPC_temporal = Constant.state_IPC_temporal;
state_power_spatial = Constant.state_power_spatial;
state_power_temporal = Constant.state_power_temporal;
state_thread_location = Constant.state_thread_location;
mig_interval = Constant.mig_interval;
mig_overhead = Constant.mig_overhead;
ignore_mig_overhead = Constant.ignore_mig_overhead;
n_step = Constant.n_step;
ignore_first_n_steps = Constant.ignore_first_n_steps;
save_profile = Constant.save_profile;   
speed_opt_method = Constant.speed_opt_method;
MDF = Constant.MDF;
use_saved_floorplan = Constant.use_saved_floorplan;
use_saved_thermal_model = Constant.use_saved_thermal_model;
constrain_n_hottest_units = Constant.constrain_n_hottest_units;
n_hottest_units = Constant.n_hottest_units;
s_tol = Constant.s_tol;
f_tol = Constant.f_tol;
thread_migration_policy = Constant.thread_migration_policy;
s_min = Constant.s_min;
n_step_per_global_opt = Constant.n_step_per_global_opt;
speed_opt_window = Constant.speed_opt_window;
use_SPEC_benchmarks = Constant.use_SPEC_benchmarks;
use_MiBench_benchmarks = Constant.use_MiBench_benchmarks;
IPC_SPEC = Constant.IPC_SPEC;
n_benchmarks = Constant.n_benchmarks;
Thread_initial_map = Constant.Thread_initial_map;
power_dir = Constant.power_dir;
save_profile_to_file = Constant.save_profile_to_file;

% Duration of each simulation time step. Must be an integer factor of 
% mig_overhead. Currently, the migration policies assume time_step = 
% mig_overhead, unless faster simulation is required, in which case, 
% time_step = mig_interval
if ignore_mig_overhead
    time_step = mig_interval;
else
    time_step = mig_overhead;
end

if strcmp(LDT_model, 'exponential')
    n_iter_LP_Exp_LDT = Constant.n_iter_LP_Exp_LDT;
    t_tol_LP_Exp_LDT = Constant.t_tol_LP_Exp_LDT;
elseif strcmp(LDT_model, 'PWL')
    n_PWL_LDT = Constant.n_PWL_LDT;
    n_iter_LP_PWL_LDT = Constant.n_iter_LP_PWL_LDT;
end


%% Directory to save floorplan and thermal model data

% Model directory is a subdirectory under './models' composed from the
% multicore configuration, i.e. number of cores in x and y directions
model_dir = sprintf('./models/%dx%d/', n_core_x, n_core_y);
% If directory does not exist, create it
if ~exist(model_dir, 'dir')
    mkdir(model_dir);
end

%% Compute floorplan and thermal model data file names

% Compose floorplan filename based on the base floorplan_name, and the
% multicore configuration.
floorplan_file_SC = ['./models/1x1/' floorplan_name '_1x1.flp'];
floorplan_file_MC = [model_dir sprintf('%s_%dx%d.flp', floorplan_name, ...
    n_core_x, n_core_y)];
% If the desired multicore floorplan does not exist, saved floorplan can't 
% be used. 
if ~exist(floorplan_file_MC, 'file')
    use_saved_floorplan = 0;
end
% Compose name of thermal model file that will store the computed thermal
% model data.
thermal_model_file = [model_dir sprintf('%s_%dx%d.mat', floorplan_name, ...
    n_core_x, n_core_y)];
% If the desired thermal model file does not exist, its saved version can't
% be used.
if ~exist(thermal_model_file, 'file')
    use_saved_thermal_model = 0;
end

%% Read single core floorplan

% The multicore floorplan is formed by creating multiple copies of a single
% core floorplan. If the single core floorplan does not exist, flag error.
if ~exist(floorplan_file_SC, 'file')
    error(['Error: Could not find single core floorplan file %s.' ...
        'Exiting\n'], floorplan_file_SC);
end
% Read the single-core floorplan and obtain geometry of each unit.
[~, W_SC, ~, ~, ~] = textread(floorplan_file_SC, '%s\t%f\t%f\t%f\t%f');
n_core = n_core_x * n_core_y;
% Number of units in single core floorplan = Number of units in each core
% of multicore flooplan
n_units_per_core = length(W_SC);
n_units_in_die = n_core * n_units_per_core;

if ~constrain_n_hottest_units
    n_hottest_units = n_units_per_core;
end

%% Create multi-core floorplan

% If saved multicore floorplan can't be used, create new floorplan.
if ~use_saved_floorplan
    create_floorplan(floorplan_file_SC, n_core_x, n_core_y, ...
        floorplan_file_MC);
end

%% Extract Hotspot model parameters

% The thermal model consists of the state-space representation of the
% thermal equivalent circuit:
%   dT/dt = A_HS T + B_HS P
%   where T and P are the temperature and power of each block (or the
%   voltage and current of each node in the equivalent circuit), and [A_HS,
%   B_HS] represents the linear system of the circuit.
% Additionally, G_die_ver, G_die_lat, G_pkg contain the thermal
% conductances of the die (vertical), die (lateral), and package (spreader,
% sink, and convection), and the corresponding capacitances are stored in
% C_die and C_pkg

% If saved thermal model can be used, load it.
if use_saved_thermal_model
   load(thermal_model_file, 'Thermal_model');
   % Copy model parameters from the saved Thermal_model structure
   A_HS = Thermal_model.A_HS; 
   B_HS = Thermal_model.B_HS; 
   G = Thermal_model.G;
   C = Thermal_model.C;
%    G_die_ver = Thermal_model.G_die_ver;
%    G_die_lat = Thermal_model.G_die_lat; 
%    C_die = Thermal_model.C_die; 
%    G_pkg = Thermal_model.G_pkg; 
%    C_pkg = Thermal_model.C_pkg;
else % If saved thermal model can't be used, create new one
    % Run Hotspot and read dumped output to obtain model parameters
%     [A_HS, B_HS, G_die_ver, G_die_lat, C_die, G_pkg, C_pkg] = ...
%         extract_Hotspot_model(floorplan_file_MC, n_units_in_die, ...
%             omit_lateral);
%     % Save the extracted model parameters to file for future use
%     Thermal_model = struct('A_HS', A_HS, 'B_HS', B_HS, 'G_die_ver', ...
%         G_die_ver, 'G_die_lat', G_die_lat, 'C_die', C_die, 'G_pkg', ...
%         G_pkg, 'C_pkg', C_pkg);
    [A_HS, B_HS, G, C] = extract_Hotspot_model(floorplan_file_MC, ...
        flat_pkg_center, omit_lateral, combine_hs_conv, r_convec, ...
        n_units_in_die);
    Thermal_model = struct('A_HS', A_HS, 'B_HS', B_HS, 'G', G, 'C', C);
    save(thermal_model_file, 'Thermal_model');
end

%% Find number of blocks

% The Hotspot thermal model contains package blocks in addition to die
% blocks. The exact number is determined by the size of the state-space
% matrices like A_HS.
[n_units_total, ~] = size(A_HS);
n_units_in_pkg = n_units_total - n_units_in_die;

% Compute identity matrix needed later for temperature computation
eye_matrix_n_units_total = eye(n_units_total);

% Indices of all blocks in a core
Block_index_core_1 = 1:n_units_per_core;
Block_index_core_1_t = Block_index_core_1';

%% Analytical thermal model for identical cores and threads
C_HS = [C.Die; C.TIM; C.Package];
G_HS = - diag(C_HS) * A_HS;

Blocks_in_die_of_core_1 = 1 : n_units_per_core;
Blocks_in_int_of_core_1 = n_units_in_die + Blocks_in_die_of_core_1;
Blocks_in_die_of_other_cores = (n_units_per_core + 1) : n_units_in_die;
Blocks_in_int_of_other_cores = ...
    (n_units_in_die + n_units_per_core + 1) : (2 * n_units_in_die);
%Blocks_in_pkg = (2 * n_units_in_die + 1) : n_units_total;
Blocks_in_pkg_except_spr = (2*n_units_in_die + 2) : n_units_total;
spr_block = 2 * n_units_in_die + 1;
Blocks_in_spr = (2*n_units_in_die + 1) : (2*n_units_in_die + 5);
Blocks_in_hs_conv = (2*n_units_in_die + 6) : n_units_total;

G_die = G_HS(Blocks_in_die_of_core_1, Blocks_in_die_of_core_1);
G_die_int = G_HS(Blocks_in_die_of_core_1, Blocks_in_int_of_core_1);
G_int = G_HS(Blocks_in_int_of_core_1, Blocks_in_int_of_core_1);
G_int_die = G_HS(Blocks_in_int_of_core_1, Blocks_in_die_of_core_1);
A_die = A_HS(Blocks_in_die_of_core_1, Blocks_in_die_of_core_1);
A_die_int = A_HS(Blocks_in_die_of_core_1, Blocks_in_int_of_core_1);
A_int = A_HS(Blocks_in_int_of_core_1, Blocks_in_int_of_core_1);
A_int_die = A_HS(Blocks_in_int_of_core_1, Blocks_in_die_of_core_1);
B_die = B_HS(Blocks_in_die_of_core_1, Blocks_in_die_of_core_1);
% B_int = B_HS(Blocks_in_int_of_core_1, Blocks_in_int_of_core_1);

for i = 1:n_units_per_core
    G_die(i, i) = G_die(i, i) + sum(G_HS(i, Blocks_in_die_of_other_cores));
    G_int(i, i) = G_int(i, i) + ...
        sum(G_HS(n_units_in_die + i, Blocks_in_int_of_other_cores));
%     A_die(i, i) = - B_die(i, i) * G_die(i, i);
%     A_int(i, i) = - B_int(i, i) * G_int(i ,i);
    A_die(i, i) = A_die(i, i) + sum(A_HS(i, Blocks_in_die_of_other_cores));
    A_int(i, i) = A_int(i, i) + ...
        sum(A_HS(n_units_in_die + i, Blocks_in_int_of_other_cores));
end

G_spr_int = G_HS(spr_block, Blocks_in_int_of_core_1);
G_int_spr = G_HS(Blocks_in_int_of_core_1, spr_block); 
g_spr = G_HS(spr_block, spr_block) - ...
    G_HS(spr_block, Blocks_in_pkg_except_spr) * ...
    (G_HS(Blocks_in_pkg_except_spr, Blocks_in_pkg_except_spr) ...
    \G_HS(Blocks_in_pkg_except_spr, spr_block));

A_spr_int = A_HS(spr_block, Blocks_in_int_of_core_1);
A_int_spr = A_HS(Blocks_in_int_of_core_1, spr_block);
a_spr = A_HS(spr_block, spr_block) - ...
    A_HS(spr_block, Blocks_in_pkg_except_spr) * ...
    (A_HS(Blocks_in_pkg_except_spr, Blocks_in_pkg_except_spr)...
    \A_HS(Blocks_in_pkg_except_spr, spr_block));

A_spr_eff = A_HS(Blocks_in_spr, Blocks_in_spr) - ...
    A_HS(Blocks_in_spr, Blocks_in_hs_conv) * ...
    (A_HS(Blocks_in_hs_conv, Blocks_in_hs_conv)...
    \A_HS(Blocks_in_hs_conv, Blocks_in_spr));

Analytical = struct(...
    'G_die', G_die, ...
    'G_die_int', G_die_int, ...
    'G_int_die', G_int_die, ...
    'G_int', G_int, ...
    'g_spr', g_spr, ...
    'G_spr_int', G_spr_int, ...
    'G_int_spr', G_int_spr, ...
    'A_die', A_die, ...
    'A_die_int', A_die_int, ...
    'A_int_die', A_int_die, ...
    'A_int', A_int, ...
    'a_spr', a_spr, ...
    'A_spr_eff', A_spr_eff, ...
    'A_spr_int', A_spr_int, ...
    'A_int_spr', A_int_spr, ...
    'B_die', B_die, ...
    'I', eye(n_units_per_core), ...
    'n_units_per_core', n_units_per_core, ...
    't_max', t_max);

%% Extract RC values for one of the cores
% R_core  = 1./G_die_ver(1:n_units_per_core, 1); 
% C_core = C_die(1:n_units_per_core);

%% Compute representative thermal time constants
% tau = R_core .* C_core;
A_HS_inv_abs_diag = diag(abs(1./A_HS));
tau_die_min = min(A_HS_inv_abs_diag(1:n_units_in_die));
tau_die_max = max(A_HS_inv_abs_diag(1:n_units_in_die));
tau_pkg = 0.35 * 140.4 * 0.5; 

%% Compute speed upper bound due to mobility degradation with temperature

% Create an inline function object that computes the speed as a function of
% maximum chip temperature by accounting for the degradation of mobility
% (and hence speed) with temperature. See read_model_parameters.m to see
% the formula used.
s_max_MDF = eval(sprintf(...
    'inline(\''((%f + %f)/(t_ini_core_max + %f))^%f\'')', ...
    t_amb_abs, t_0_C_in_K, t_0_C_in_K, MDF));



%--------------------------------------------------------------------------
%% Read power data
%--------------------------------------------------------------------------
% If the power data file specified in the Constant strucuture does not
% exist, flag error and quit.
if ~exist(power_data_file, 'file')
    error('Error: Could not find power data file %s. Exiting\n', ...
        power_data_file);
end
% If file exists, extract normalized dynamic, active leakage, and standby
% leakage power values from file. To see the file format for the
% power_data_file, type 'help read_power_data'.
if use_SPEC_benchmarks || use_MiBench_benchmarks
    [P_d_norm, P_a_norm, P_s_norm, T] = read_power_data(...
        power_dir, n_units_per_core, 0, 1);
else
    [P_d_norm, P_a_norm, P_s_norm, T] = read_power_data(...
        power_data_file, n_units_per_core);
end
if use_SPEC_benchmarks
    P_s_norm = mean(P_s_norm, 3);
end

if Constant.use_MiBench_benchmarks
        power_data_str = {'MiBench/basicmath_large.txt'; ...
        'MiBench/qsort_large.txt'; ...
        'MiBench/ffti_large.txt'; ...
        'MiBench/gsm_decode_large.txt'};
    
    P_d_ini = load(power_data_str{1});
    max_inst = Constant.max_inst;
    power_data_len = length(power_data_str);
    P_d_norm = zeros(n_units_per_core, power_data_len, ...
        max(Constant.max_inst, size(P_d_ini,1)));
    
    for c = 1:power_data_len
        P_d_ini = load(power_data_str{c});
        P_d_ini_len = size(P_d_ini,1);
        if P_d_ini_len >= max_inst
            P_d_norm(:,c,:) = P_d_ini(1:max_inst);
        else
            tot_len = 0;
            while true
                if (max_inst - tot_len) > P_d_ini_len
                    P_d_norm(:,c,tot_len + 1:tot_len + P_d_ini_len) ...
                        = P_d_ini(:,:)';
                    tot_len = tot_len + P_d_ini_len;
                else
                    P_d_norm(:,c,tot_len + 1:end) ...
                        = P_d_ini(1:(max_inst - tot_len),:)';
                    break;
                end
            end
        end
    end  
end
    P_d_norm = P_d_norm/max(sum(mean(P_d_norm,3),1));

T = T - t_amb_abs;
% Make sure that the value of maximum temperature used in file is same as
% that specified in the Constant strucuture.
if T(end) ~= t_max
    error(['Error: Inconsistency between maximum temperature ' ...
        'chosen and corresponding value read from power data file\n']);
end
if T(1) ~= 0
    error(['Error: Inconsistency between ambient temperature ' ...
        'chosen and corresponding value read from power data file\n']);
end

%% Scale normalized power data to get (mean) per core power data
P_a_mean_LDT = p_a_at_Tmax_for_SC * P_a_norm / n_core;
P_s_mean_LDT = p_s_at_Tmax_for_SC * P_s_norm / n_core;
if strcmp(LDT_model, 'constant') 
    if use_SPEC_benchmarks
        P_a_norm = squeeze(P_a_norm(:, end, :));
    else
        P_a_norm = P_a_norm(:, end);
        P_s_norm = P_s_norm(:, end);
    end
end
P_d_mean = p_d_at_fmax_for_SC * P_d_norm / n_core;
P_a_mean = p_a_at_Tmax_for_SC * P_a_norm / n_core;
P_s_mean = p_s_at_Tmax_for_SC * P_s_norm / n_core;

%% Perform piece-wise linearization of standby LDT

% The exponential power model also computes the linear model coefficients,
% which are used to compute the starting solution for faster convergence.
% Compute standby leakage power coefficients by finding a linear upper
% bound. Assuming P(T) is convex, choose slope so that linear aprprox.
% lies above convex curve
% P = C1 + C2 * T
    
P_s_coeff_lin_LDT = zeros(n_units_per_core, 2);
rmse = zeros(n_units_per_core, 1);
for i_unit = 1:n_units_per_core         
    C1 = P_s_mean_LDT(i_unit, 1);
    C2 = (P_s_mean_LDT(i_unit, end) - P_s_mean_LDT(i_unit, 1))...
        /(T(end) - T(1));
    P_s_coeff_lin_LDT(i_unit, :) = [C1 C2];
    P_s_mean_pred = C1 + C2 * T;
    rmse(i_unit) = mean(abs(P_s_mean_LDT(i_unit, :)' - P_s_mean_pred));
end    
if strcmp(LDT_model, 'linear')
    P_s_coeff_LDT = P_s_coeff_lin_LDT;
end

if strcmp(LDT_model, 'PWL')
    P_s_coeff_LDT = zeros(n_units_per_core, 2, n_PWL_LDT);
    % Split temperature range into multiple regions
    [I_split, T_split] = split_vector(T, n_PWL_LDT);
    T_split = [T_split; T(end)];
    % Create a matrix with multiple copies of T_split. This will be used
    % later to find the temperature region to which each die block belongs
    % to by comparing a repeated form of the temperature vector to this
    % matrix. See temperature_step_response_PWL_LDT.m for more details.
    T_split_1 = T_split';
    T_split_1(end) = T_split(end) + 5; % Hack. 
    T_split_rep = repmat(T_split_1, n_units_in_die, 1);
    % Compute standby power coefficients by doing a PWL fit
    for i_unit = 1:n_units_per_core
        [C1, C2] = piecewise_linear_fit(...
            T, P_s_mean(i_unit,:)', I_split, n_PWL_LDT, 0);
        P_s_coeff_LDT(i_unit, 1, :) = C1;
        P_s_coeff_LDT(i_unit, 2, :) = C2;
    end    
end

% Compute standby leakage power coefficients by doing an exp. fit
% P = C1 + C2 * exp(C3 * T). This is to be done regardless of the LDT model
% chosen because the temperature computation is done with the most accurate
% LDT model available. The selected LDT_model is used only in finding the
% optimum speed of each core.
P_s_coeff_exp_LDT = zeros(n_units_per_core, 3);
coeff0 = zeros(3, 1);
exp_LDT_model_fun = ...
    @(coeff, T1) coeff(1) + coeff(2) * exp(coeff(3) * T1);
rmse = zeros(n_units_per_core, 1);
for i_unit = 1:n_units_per_core
    % Initial guess for model coefficients
    coeff0(1) = P_s_mean_LDT(i_unit, 1)';
    coeff0(2) = coeff0(1);
    coeff0(3) = ...
        log((P_s_mean_LDT(i_unit, end) - coeff0(1))/coeff0(2)) / T(end);
    coeff = ...
        nlinfit(T, P_s_mean_LDT(i_unit, :)', exp_LDT_model_fun, coeff0);
    P_s_coeff_exp_LDT(i_unit, :) = coeff';
    P_s_mean_pred = coeff(1) + coeff(2) * exp(coeff(3) * T);
    rmse(i_unit) = mean(abs(P_s_mean_LDT(i_unit, :)' - P_s_mean_pred));
end
if strcmp(LDT_model, 'exponential')
    P_s_coeff_LDT = P_s_coeff_exp_LDT;
end

%% Generate spatial thread IPC profiles
% Set state of random number of generator to specified value.
rand(rand_gen_method, state_IPC_spatial);
if use_SPEC_benchmarks
    IPC = IPC_SPEC(Thread_initial_map);
else
    % Generate n_thread random numbers uniformly distributed over 
    % [1-IMR_IPC_spatial/2, 1+IMR_IPC_spatial/2]. These constitute
    % (multiplicative) deviations about the IPC_mean.
    IPC_dev = 1 - IMR_IPC_spatial/2 + IMR_IPC_spatial*rand(n_thread, 1);
    % Use these values to get the IPCs for each thread
    IPC = IPC_dev * IPC_mean;
end

%% Generate spatial dynamic power profiles
rand(rand_gen_method, state_power_spatial);
if use_SPEC_benchmarks || use_MiBench_benchmarks
    P_d = P_d_mean(:, Thread_initial_map,:);
else
    P_d_dev = 1 - IMR_power_spatial/2 ...
        + IMR_power_spatial*rand(n_units_per_core, n_thread);
    P_d = P_d_dev .* repmat(P_d_mean, 1, n_thread);
end

%% Generate spatial active leakage power profiles
if use_SPEC_benchmarks
    if strcmp(LDT_model, 'constant')
        P_a = P_a_mean(:, Thread_initial_map);
    else
        P_a = P_a_mean(:, :, Thread_initial_map);
    end
else
    P_a_dev = 1 - IMR_power_spatial/2 ...
        + IMR_power_spatial*rand(n_units_per_core, n_thread);
    if strcmp(LDT_model, 'constant')
        P_a = P_a_dev .* repmat(P_a_mean, 1, n_thread);
    end
end
% Compute standby leakage power coefficients by finding a linear upper
% bound. Assuming P(T) is convex, choose slope so that linear 
% approx. lies above convex curve
% P = C1 + C2 * T
% The exponential power model also computes the linear model coeff.'s,
% which are used to find the starting solution for faster convergence.
P_a_coeff_lin_LDT = zeros(n_units_per_core, n_thread, 2);
rmse = zeros(n_thread, n_units_per_core);
for i_thread = 1:n_thread
    for i_unit = 1:n_units_per_core
        if use_SPEC_benchmarks
            %P_a_lin = squeeze(P_a_mean_LDT(i_unit, :, i_thread))';
            P_a_lin = squeeze(P_a(i_unit, :, i_thread))';
        else
            P_a_lin = (P_a_mean_LDT(i_unit, :) * ...
                P_a_dev(i_unit, i_thread))';
        end
        C1 = P_a_lin(1);
        C2 = (P_a_lin(end) - P_a_lin(1)) / (T(end) - T(1));
        P_a_coeff_lin_LDT(i_unit, i_thread, :) = [C1 C2];
        P_a_mean_pred = C1 + C2 * T;
        rmse(i_thread, i_unit) = mean(abs(P_a_lin - P_a_mean_pred));
    end
end
if strcmp(LDT_model, 'linear')
    P_a_coeff_LDT = P_a_coeff_lin_LDT;
end

if strcmp(LDT_model, 'PWL')
    P_a_coeff_LDT = zeros(n_units_per_core, n_thread, 2, n_PWL_LDT);
    for i_thread = 1:n_thread
        for i_unit = 1:n_units_per_core
            % Add code for handling SPEC benchmarks case
            if use_SPEC_benchmarks
                P_a_PWL = squeeze(P_a(i_unit, :, i_thread))';
            else
                P_a_PWL = (P_a_mean(i_unit, :) * ...
                    P_a_dev(i_unit, i_thread))';
            end
            [C1, C2] = piecewise_linear_fit(...
                T, P_a_PWL, I_split, n_PWL_LDT, 0);
            P_a_coeff_LDT(i_unit, i_thread, 1, :) = C1;
            P_a_coeff_LDT(i_unit, i_thread, 2, :) = C2;
        end
    end
end

% Compute standby leakage power coefficients by doing an exp. fit
% P = C1 + C2 * exp(C3 * T). This is to be done regardless of the LDT model
% chosen because the temperature computation is done with the most accurate
% LDT model available. The selected LDT_model is used only in finding the
% optimum speed of each core.
P_a_coeff_exp_LDT = zeros(n_units_per_core, n_thread, 3);
coeff0 = zeros(3, 1);
exp_LDT_model_fun = ...
    @(coeff, T) coeff(1) + coeff(2) * exp(coeff(3) * T);
rmse = zeros(n_thread, n_units_per_core);
for i_thread = 1:n_thread
    for i_unit = 1:n_units_per_core
        if use_SPEC_benchmarks
            P_a_exp = squeeze(P_a(i_unit, :, i_thread))';
        else
            P_a_exp = (P_a_mean_LDT(i_unit, :) * ...
                P_a_dev(i_unit, i_thread))';
        end
        % Initial guess for model coefficients
        coeff0(1) = P_a_exp(1);
        coeff0(2) = coeff0(1);
        coeff0(3) = ...
            log((P_a_exp(end) - coeff0(1))/coeff0(2)) / T(end);
        coeff = nlinfit(T, P_a_exp, exp_LDT_model_fun, coeff0);
        P_a_coeff_exp_LDT(i_unit, i_thread, :) = coeff';
        P_a_mean_pred = coeff(1) + coeff(2) * exp(coeff(3) * T);
        rmse(i_thread, i_unit) = mean(abs(P_a_exp - P_a_mean_pred));
    end
end    
if strcmp(LDT_model, 'exponential')
    P_a_coeff_LDT = P_a_coeff_exp_LDT;
end

%% Set number of threads that are rotated between successive thread 
% migration intervals. 

% See sim_rotation_mig_policy.m to see where this is used.
n_max_thread_per_mig = n_thread;

%% Temporal power variations are computed as random perturbations about the
% mean power computed above. These are done every
% n_step_per_power_variation time steps, which is chosen so that the time
% between successive power changes is of the order of 0.1 times the
% smallest die thermal time constant.
n_step_per_power_variation = max(1, floor(0.1 * tau_die_min / time_step));

%% Create a one's vector that will be useful later in computing linear
% indices for the leakage power coefficient matrices
Ones = ones(n_units_per_core, 1);

%% Initial mapping of threads to cores

% An array that stores the core row (Y) index for each thread 
Row = zeros(n_thread, 1);
% An array that stores the core column (X) index for each thread 
Col = zeros(n_thread, 1);
% Thread index
cur_thread = 1;
% Matrix of same size as multicore processor. It stores the thread mapping,
% i.e. the value of On_core(r, c) is the id of the thread that is mapped to
% core row r and column c, with an id of 0 indicating an empty core.
On_core = zeros(n_core_y, n_core_x);
% Set state of random number generator
rand(rand_gen_method, state_thread_location);
% For each thread, randomly generate core location, and assign it
if use_SPEC_benchmarks
    row = 1;
    col = 1;
end
while cur_thread <= n_thread
    if ~use_SPEC_benchmarks
        % row ~ U[1, n_core_y], col ~ U[1, n_core_x]; row, col are integers
        row = 1 + floor(n_core_y * 0.99 * rand(1, 1));
        col = 1 + floor(n_core_x *0.99 * rand(1, 1));
    end
    Row(cur_thread) = row;
    Col(cur_thread) = col;   
    % Check if this random combination occured before
    is_duplicate = 0;
    if cur_thread > 1
        % For all threads seen before
        for prev_thread = 1:(cur_thread-1)
            % If prev_thread was assigned same location as cur_thread
            if (row == Row(prev_thread)) && (col == Col(prev_thread))
                % Mark this location as duplicate
                is_duplicate = 1;
                break;
            end
        end
    end
            
    if ~is_duplicate % If not duplicate, store location in On_core matrix
        On_core(row, col) = cur_thread;
        cur_thread = cur_thread + 1;
    end
    if use_SPEC_benchmarks
        col = col + 1;
        if col > n_core_x
            col = 1;
            row = row + 1;
        end
    end    
end

%% Power vector/coeffs and IPC's corresponding to the above initial 
% thread mapping.  
Block = 1:n_units_per_core;
P_d_ini = zeros(n_units_total, max_inst);
switch LDT_model
    case 'constant',
        P_l_ini = zeros(n_units_total, 1);
    case 'linear',
        P_l_coeff_LDT_ini = zeros(n_units_total, 2);
    case 'PWL',
        P_l_coeff_LDT_ini = zeros(n_units_total, 2, n_PWL_LDT);
    case 'exponential',
        P_l_coeff_LDT_ini = zeros(n_units_total, 3);
end
P_l_coeff_lin_LDT_ini = zeros(n_units_total, 2);
P_l_coeff_exp_LDT_ini = zeros(n_units_total, 3);

% Zeta_ini = zeros(n_units_in_die, 1);
% P_d_prime_ini = zeros(n_units_in_die, 1);
% P_l_prime_ini = zeros(n_units_in_die, 1);
% G = 0;

IPC_ini = zeros(n_core, 1);

core = 1;
for core_y = 1:n_core_y
    for core_x = 1:n_core_x
        thread = On_core(core_y, core_x);
        if thread ~= 0
            IPC_ini(core) = IPC(thread);
            P_d_ini(Block,:) = P_d(:, thread,:);
            P_l_coeff_exp_LDT_ini(Block, :) = ...
                squeeze(P_a_coeff_exp_LDT(:, thread, :));                      
            P_l_coeff_lin_LDT_ini(Block, :) = ...
                squeeze(P_a_coeff_lin_LDT(:, thread, :));
                        
            switch LDT_model
                case 'constant',
                    P_l_ini(Block) = P_a(:, thread);
                case 'linear', 
                    P_l_coeff_LDT_ini(Block, :) = ...
                        squeeze(P_a_coeff_LDT(:, thread, :));
                case 'PWL',
                    P_l_coeff_LDT_ini(Block, :, :) = ...
                        squeeze(P_a_coeff_LDT(:, thread, :, :));
                case 'exponential',
                    P_l_coeff_LDT_ini(Block, :) = ...
                        squeeze(P_a_coeff_LDT(:, thread, :));                    
            end
        else
            P_d_ini(Block) = 0;
            P_l_coeff_exp_LDT_ini(Block, :) = P_s_coeff_exp_LDT;
            P_l_coeff_lin_LDT_ini(Block, :) = P_s_coeff_lin_LDT;

            switch LDT_model
                case 'constant',
                    P_l_ini(Block) = P_s_mean;
                case 'linear',
                    P_l_coeff_LDT_ini(Block, :) = P_s_coeff_LDT;
                case 'PWL',
                    P_l_coeff_LDT_ini(Block, :, :) = P_s_coeff_LDT;
                case 'exponential',
                    P_l_coeff_LDT_ini(Block, :) = P_s_coeff_LDT;                    
            end
        end
        % Analytical power-thermal model parameters
%         Zeta_ini(Block) = 1./ ... 
%             (1 - P_l_coeff_lin_LDT_ini(Block, 2) .* R_die_int_ver(Block));
%         P_d_prime_ini(Block) = Zeta_ini(Block) .* P_d_ini(Block);
%         P_l_prime_ini(Block) = Zeta_ini(Block) .* ...
%             P_l_coeff_lin_LDT_ini(Block, 1);
%         G = G + sum((Zeta_ini(Block) - 1) ./ R_die_int_ver(Block));
        
        Block = Block + n_units_per_core;
        core = core + 1;
    end
end
% Note that if LDT is accounted for, then we need to first know the initial
% temperature to compute the static power for the initial thread mapping.

%% Compute approximate thermal model
N_hs_conv = 2*n_units_in_die + 1 + 4 + (1:9);
C_all = [C.Die; C.TIM; C.Package];
Hotspot_Approx = compute_approx_thermal_model(A_HS, B_HS, C_all, ...
    N_hs_conv);

%% Precomputed matrices to be used later in temperature computation
% The power vector P (one element for each block in die+package system) as 
% a function of the speed vector S (one element for each core) is:
%   P(S) = P_d(S) + P_a + P_s, where P_d, P_a, and P_s, are the dynamic,
% active leakage, and the standby leakage components, respectively. The
% dynamic power itself can be expressed as
%   P_d = diag(P_d_max) * D * S
% where D is a binary matrix chosen such that for each unit k of the chip
% that lies in core i, P_d(k) = P_d_max(k) * S(i), for each unit j of the
% chip that lies in the pacakge, P_d(j) = 0.
%
% The temperature at the end of a time step T in response to the step power
% vector input P at the beginning of the time step, for a given initial
% temperature T_ini at the beginning of the time step, can be computed from
% state-space concepts as follows:
%   T(S) = expm(A_HS * time_step_duration) * T_ini + inv(A_HS) *
%       (expm(A_HS) - I) * B_HS * P(S)
% where S is the speed vector consisting of speed of all cores, and P(S)
% and T(S) are the corresponding power and temperature vectors, and expm
% is the matrix exponential operation. 
%
% This can be written as
%   T(S) = M * T_ini + R * diag(P_d_max) * D + R * (P_a + P_s), where
%       M = expm(A_HS * time_step_duration),
%       R = inv(A_HS) * (expm(A_HS) * time_step_duration - I) * B_HS
% As the above computation is repeated at least 10^5 times during transient
% simulations, it is worthwhile to precompute, M, R, and D.

Precomputed = struct(...
    'M', zeros(n_units_total, n_units_total), ...
    'M_glob', zeros(n_units_total, n_units_total), ...
    'M_lin', zeros(n_units_total, n_units_total), ...
    'R_lin', zeros(n_units_total, n_units_total), ...
    'R', zeros(n_units_total, n_units_total), ... 
    'R_glob', zeros(n_units_total, n_units_total), ...     
    'D', zeros(n_units_total, n_core), ...
    'I_n_hottest', zeros(n_hottest_units * n_core, 1), ...
    'R_star_P_s', zeros(n_units_total, 1), ...
    'T_max_minus_R_star_P_s', zeros(n_units_total, 1), ...
    'R_star_P_s_lin', zeros(n_units_total, 1), ...
    'T_max_minus_R_star_P_s_lin', zeros(n_units_total, 1), ...
    'P_d_star_D', zeros(n_units_total, n_core), ...
    'R_star_P_d_star_D', zeros(n_units_total, n_core));

if ~constrain_n_hottest_units
    % All die units
    Precomputed.I_n_hottest = 1:n_units_in_die;
end

switch LDT_model
    case 'constant',
        Precomputed.M = expm(A_HS * time_step);
        Precomputed.R = (A_HS\ ...
            (Precomputed.M - eye(n_units_total))) * B_HS;
        if strcmp(speed_opt_window, 'local_global') || ...
                strcmp(speed_opt_window, 'global')
            Precomputed.M_global = ...
                expm(A_HS * time_step * n_step_per_global_opt);
            Precomputed.R_global = (A_HS\ ...
                (Precomputed.M_global - eye(n_units_total))) * B_HS;    
        end
    case'linear',
        A_HS_prime = A_HS + B_HS * diag(P_l_coeff_LDT_ini(:,2));            
        Precomputed.M = expm(A_HS_prime * time_step);
        Precomputed.R = (A_HS_prime\ ...
            (Precomputed.M - eye(n_units_total)))* B_HS;
        if strcmp(speed_opt_window, 'local_global') || ...
                strcmp(speed_opt_window, 'global')
            Precomputed.M_global = ...
                expm(A_HS_prime * time_step * n_step_per_global_opt);
            Precomputed.R_global = (A_HS_prime\ ...
                (Precomputed.M_global - eye(n_units_total))) * B_HS;   
        end
end
A_HS_prime = A_HS + B_HS * diag(P_l_coeff_lin_LDT_ini(:,2));            
Precomputed.M_lin = expm(A_HS_prime * time_step);
Precomputed.R_lin = A_HS_prime\ ...
    (Precomputed.M_lin - eye(n_units_total)) * B_HS;

Precomputed.D = zeros(n_units_total, n_core);
Row_Index = 1:n_units_per_core;
for core = 1:n_core
    D_row = [zeros(1, (core-1)), 1, zeros(1, n_core-core)];
    Precomputed.D(Row_Index, :) = repmat(D_row, n_units_per_core, 1);
    Row_Index = Row_Index + n_units_per_core;
end

switch LDT_model 
    case 'constant',
        switch speed_opt_method
            case 'bin_search',    
                Precomputed.R_star_P_s = Precomputed.R * P_l_ini;                             
            case 'LP_NLP',
                Precomputed.T_max_minus_R_star_P_s = - Precomputed.R * ...
                    P_l_ini + t_max * ones(n_units_total, 1);
        end
    case 'linear',
        switch speed_opt_method
            case 'bin_search',    
                Precomputed.R_star_P_s = Precomputed.R * ...
                    P_l_coeff_LDT_ini(:,1);           
            case 'LP_NLP',
                Precomputed.T_max_minus_R_star_P_s = - Precomputed.R * ...
                    P_l_coeff_LDT_ini(:,1) + ...
                    t_max * ones(n_units_total, 1);
        end        
end

switch speed_opt_method
    case 'bin_search',    
        Precomputed.R_star_P_s_lin = Precomputed.R_lin * ...
            P_l_coeff_lin_LDT_ini(:,1);     
    case 'LP_NLP',
        Precomputed.T_max_minus_R_star_P_s_lin = - Precomputed.R_lin * ...
            P_l_coeff_lin_LDT_ini(:,1) + ...
            t_max * ones(n_units_total, 1);
end

%% Precompute some matrices needed for temperature computation using the 
% approximate thermal model
n_die_int_spr = N_hs_conv(1)-1; 
N_die_int_spr = 1:n_die_int_spr;
A1_prime = Hotspot_Approx.A1 + Hotspot_Approx.B1 * ...
    diag(P_l_coeff_lin_LDT_ini(N_die_int_spr, 2));
tau_conv_eff = 1/(Hotspot_Approx.A3 * (A1_prime\ ...
    Hotspot_Approx.A2) - Hotspot_Approx.a4);
inv_A1 = inv(A1_prime);
W1 = -inv_A1 * Hotspot_Approx.B1;
W2 = -inv_A1 * Hotspot_Approx.A2;
W3 = W2 * Hotspot_Approx.A3 * W1;
W4 = -Hotspot_Approx.A3 * W1;
Blocks_in_die = 1:n_units_in_die;
Precomputed.Approx_Model = struct(...
    'W1', W1(Blocks_in_die, Blocks_in_die), ...
    'W2', W2(Blocks_in_die), ...
    'W3', W3(Blocks_in_die, Blocks_in_die), ...
    'W4', W4(Blocks_in_die), ...
    'tau_conv_eff', tau_conv_eff);

% To compute die temperature as a function of speed

%% Compute some options needed for LP/NLP speed optimization method
if strcmp(speed_opt_method, 'LP_NLP')
    % Speed lower bound for each core
    LB = s_min * ones(n_core, 1);
    % Speed upper bound for each core
    UB = ones(n_core, 1);
    % Options for linear program solver
    options = optimset('Display', 'off');
    
    LP_NLP = struct(...
    'LB', LB, ...
    'UB', UB, ...
    'options', options);
end

%% Create data structures
Floorplan = struct(...
    'name', floorplan_name, ...
    'single_core_FP_file', floorplan_file_SC, ...
    'multi_core_FP_file', floorplan_file_MC);  

Die = struct(...
    'n_core_x', n_core_x, ...
    'n_core_y', n_core_y, ...
    'n_core', n_core, ...
    'omit_lateral', omit_lateral, ...
    'Floorplan', Floorplan, ...
    's_min', s_min);

N_units = struct(...
    'per_core', n_units_per_core, ...
    'in_die', n_units_in_die, ...
    'in_pkg', n_units_in_pkg, ...
    'total', n_units_total);

Hotspot = struct(...
    'A', A_HS, ...
    'B', B_HS, ...
    'G', G_HS, ...
    'C', C_HS, ...    
    'Approx', Hotspot_Approx);

Time_constant_die = struct(...
    'min', tau_die_min, ...
    'max', tau_die_max);

Time_constant = struct(...
    'Die', Time_constant_die, ...
    'pkg', tau_pkg);

Thermal = struct(...
    'N_units', N_units, ...
    'I', eye_matrix_n_units_total, ...
    'Hotspot', Hotspot, ...
    't_max_abs', t_max_abs, ...
    't_amb_abs', t_amb_abs, ...
    't_max', t_max, ...
    't_0_C_in_K', t_0_C_in_K, ...
    's_max_MDF', s_max_MDF, ...
    'MDF', MDF, ...
    'Time_constant', Time_constant, ...
    'Precomputed', Precomputed, ...
    'Block_index_core_1', Block_index_core_1, ...
    'Block_index_core_1_t', Block_index_core_1_t, ...
    'G', G, ...
    'C', C, ...
    'G_HS', G_HS, ...
    'Analytical', Analytical);

Mean_power = struct(...
    'Dynamic', P_d_mean, ...
    'ActLeak', P_a_mean, ...
    'StdLeak', P_s_mean);

Mean_power.StdLeakCoeffExp = P_s_coeff_exp_LDT;
Mean_power.StdLeakCoeffLin = P_s_coeff_lin_LDT;
if ~strcmp(LDT_model, 'constant')
    Mean_power.StdLeakCoeff = P_s_coeff_LDT;
end
    
Power = struct(...
    'LDT_model', LDT_model, ...
    'Temperature_range', T, ...
    'power_data_file', power_data_file, ...
    'Mean', Mean_power);

if strcmp(LDT_model, 'PWL')
    Power.n_PWL_LDT = n_PWL_LDT;
    Power.T_split_rep = T_split_rep;
    Power.T_split = T_split_1;
end
    
IMR = struct(...
     'IPC_spatial', IMR_IPC_spatial, ...
     'IPC_temporal', IMR_IPC_temporal, ...
     'power_spatial', IMR_power_spatial, ...
     'power_temporal', IMR_power_temporal);
 
State = struct(...
     'IPC_spatial', state_IPC_spatial, ...
     'IPC_temporal', state_IPC_temporal, ...
     'power_spatial', state_power_spatial, ...
     'power_temporal', state_power_temporal);
 
Location = struct(...
    'Row', Row, ...
    'Col', Col);

switch LDT_model
    case 'constant',
        P_ini = struct(...
            'Dynamic', P_d_ini, ...
            'Leakage', P_l_ini);
        Thread_power_profiles = struct(...
            'Dynamic', P_d, ...
            'ActLeak', P_a);
    case {'linear', 'PWL', 'exponential'}
        P_ini = struct(...
            'Dynamic', P_d_ini, ...
            'LeakageCoeff', P_l_coeff_LDT_ini);
        Thread_power_profiles = struct(...
            'Dynamic', P_d, ...
            'ActLeakCoeff', P_a_coeff_LDT);
end
P_ini.LeakageCoeffLin = P_l_coeff_lin_LDT_ini;
P_ini.LeakageCoeffExp = P_l_coeff_exp_LDT_ini;
Thread_power_profiles.ActLeakCoeffLin = P_a_coeff_lin_LDT;
Thread_power_profiles.ActLeakCoeffExp = P_a_coeff_exp_LDT;

% Analytical = struct(...
%     'Zeta_ini', Zeta_ini, ...
%     'P_d_prime_ini', P_d_prime_ini, ...
%     'P_l_prime_ini', P_l_prime_ini, ...
%     'G', G, ...
%     'R_p', sum(R_ver.Pkg), ...
%     'R_die_int', R_die_int_ver);

Initial_map = struct(...
    'On_core', On_core, ...
    'Location', Location, ...
    'Power', P_ini, ...
    'IPC', IPC_ini);

%..........................................................................
% Thread = struct(...
%     'IMR', IMR, ...
%     'rand_gen_method', rand_gen_method, ...
%     'State', State, ...
%     'n_thread', n_thread, ...
%     'IPC_mean', IPC_mean, ...
%     'IPC', IPC, ...
%     'Initial_map', Initial_map, ...
%     'Power', Thread_power_profiles, ...
%     'Ones', Ones);
% 


Thread.IMR= IMR, ...
Thread.rand_gen_method= rand_gen_method, ...
Thread.State= State, ...
Thread.n_thread= n_thread, ...
Thread.IPC_mean= IPC_mean, ...
Thread.IPC= IPC, ...
Thread.Initial_map= Initial_map, ...
Thread.Power= Thread_power_profiles, ...
Thread.Ones= Ones
%..........................................................................

% Simulation = struct(...
%     'n_max_thread_per_mig', n_max_thread_per_mig, ...
%     'mig_interval', mig_interval, ...
%     'mig_overhead', mig_overhead, ...
%     'ignore_mig_overhead', ignore_mig_overhead, ...
%     'n_step', n_step, ...
%     'time_step', time_step, ...
%     'save_profile', save_profile, ...
%     'save_profile_to_file', save_profile_to_file, ...
%     'ignore_first_n_steps', ignore_first_n_steps, ...
%     'speed_opt_method', speed_opt_method, ...
%     'n_step_per_power_variation', n_step_per_power_variation, ...
%     'constrain_n_hottest_units', constrain_n_hottest_units, ...
%     'n_hottest_units', n_hottest_units, ...
%     'n_step_per_global_opt', n_step_per_global_opt, ...
%     'speed_opt_window', speed_opt_window);

% maximum number of thereads per migrations 
Simulation.n_max_thread_per_mig= n_max_thread_per_mig; 
% Migartion interval in Sec ??
Simulation.mig_interval= mig_interval
% Migration overhead specified in 
Simulation.mig_overhead= mig_overhead
% Ignore migartion overhead when this flag is set (1)
Simulation.ignore_mig_overhead= ignore_mig_overhead
% Number of steps temperature to be computed 
Simulation.n_step= n_step
% Sampling time of each step specified in ??
Simulation.time_step= time_step

% Save profile of temperature when thsi flag is set
Simulation.save_profile= save_profile
% Save profile to a file when thsi si set (what is the diff with teh
% above?)
Simulation.save_profile_to_file= save_profile_to_file

%Ignoe first n steps, specified in time_step
Simulation.ignore_first_n_steps= ignore_first_n_steps
%Speed optimization method specied using string 
Simulation.speed_opt_method= speed_opt_method
% Number of steps per power variations specified in int??
Simulation.n_step_per_power_variation= n_step_per_power_variation
% Number of units that is contrained for hotspot, specified in int
Simulation.constrain_n_hottest_units= constrain_n_hottest_units
% Number of hotest unitm, specified in int
Simulation.n_hottest_units= n_hottest_units
% Number of steps used per global optimization
Simulation.n_step_per_global_opt= n_step_per_global_opt
% Speed optimization window, specified in int
Simulation.speed_opt_window= speed_opt_window

%..........................................................................

if strcmp(LDT_model, 'exponential') 
    LP_Exp_LDT = struct(...
        'n_iter', n_iter_LP_Exp_LDT, ...
        't_tol', t_tol_LP_Exp_LDT);
    Simulation.LP_Exp_LDT = LP_Exp_LDT;
elseif strcmp(LDT_model, 'PWL')
    Simulation.n_iter_LP_PWL_LDT = n_iter_LP_PWL_LDT;
end

Throttling = [];
if strcmp(speed_opt_method, 'bin_search')
    Binary_search = struct(...
        's_tol', s_tol, ...
        'f_tol', f_tol);
    Throttling = struct('Binary_search', Binary_search);
elseif strcmp(speed_opt_method, 'LP_NLP');
    Throttling = struct('LP_NLP', LP_NLP);
end

Policy = struct(...
    'Throttling', Throttling, ...
    'thread_migration', thread_migration_policy);

%--------------------------------------------------------------------------
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================

% Script for testing temperature_step_response_PWL_LDT

fprintf('Computing system models....\t');
[Sim_params, Power1]= collect_sim_params(Params)
fprintf('Finished\n\n');



%     [Die, Power, Thermal, Thread, Simulation, Policy]= ...
%     compute_model_parameters(Constant)
   
% %% Compute parameters
% 
%     n_core = Params.n_core_x*Params.n_core_y;
%     core_units = Thermal.N_units.per_core;
%     total_units = Thermal.N_units.total;
%     vt = 0.3/vdd_max;
% 
%     vdd_range = 0:0.1:1;
%     vdd_len = length(vdd_range);
% 
%     P_s0 = Power.Mean.ActLeak(:,:,1)/n_core;
%     
%     % Smoothen leakage power w.r.t temperature
%     for i = 1:core_units
%         P_s0(i,:) = smooth(P_s0(i,:));
%     end
%     
%     temp_len = length(Power.Temperature_range);
%     P_s = zeros(n_core*core_units, temp_len, vdd_len);
%     
%     stream = RandStream('mrg32k3a','seed', 4839);
%                 
%     max_inst = Constant.max_inst;
%     P_d = zeros(core_units, n_core, max_inst);
%     if strcmp(dtm_op,'scaling') && n_task <= 4
%          P_d = Thread.Power.Dynamic/n_core;
%      else
%         for task = 1:n_task
%             P_d(:,task,:) = squeeze(Thread.Power.Dynamic(:,1,:)) ...
%                 * (0.25 + 1.0*rand(stream));
%         end
%      end
%     
%     for core = 1:n_core
%         die_units = (core-1)*core_units+1:core*core_units;
%         P_s0_rand = ...
%             P_s0.*repmat((0.85 + 0.3*rand(stream,core_units,1)),1,temp_len);
%         vdd_index = 0;
%         for vdd = vdd_range
%             vdd_index = vdd_index + 1;
%             P_s(die_units,:,vdd_index) = P_s0_rand*exp(vdd);
%         end
%     end
%     
% %% Build speed-voltage-temperature lookup table    
%     K = 273.15; 
%     spd_range = 0:0.1:1;
%     vdd_tab = zeros(length(spd_range),temp_len);
%     
%     spd_ind = 0;
%     for curr_spd = spd_range
%         spd_ind = spd_ind + 1;
%         temp_ind = 0;
%         for curr_temp = Power.Temperature_range' + Thermal.t_amb_abs + K
%             temp_ind = temp_ind + 1;
%             vdd_tab(spd_ind,temp_ind) = ...
%                 real(fsolve(@error_spd, 1, optimset('Display','off')));
%         end
%     end
%     
% %%    
%     Migration = struct();
%     
%     if strcmp(dtm_op, 'migration')
%         G_app_inv = zeros(n_core, core_units, core_units);
%         G_coeff_Tspr = zeros(n_core, core_units);
%         G_ldt = zeros(core,core_units);
%         
%         for core = 1:n_core
%             G_die = Thermal.Analytical.G_die;
%             G_int =  Thermal.Analytical.G_int;
%             G_ldt(n_core,:) = P_s0(:,5) ...
%                 .*(0.85 + 0.3*rand(stream,core_units,1));
%             G_die_int =  Thermal.Analytical.G_int_die;
%             G_int_spr =  Thermal.Analytical.G_int_spr;
%             
%             G_app_inv(core,:,:) = inv(G_die - ...
%                 G_die_int*inv(G_int)*G_die_int - diag(G_ldt(n_core,:)));
%             G_coeff_Tspr(core,:,:) = G_die_int*inv(G_int)*G_int_spr;
%         end
%         
%         time_const = Thermal.Time_constant.pkg;
%         dec_rate = exp(-t_step/time_const);
%         inc_rate = 1-dec_rate;
%         
%         Migration.G_app_inv = G_app_inv;
%         Migration.G_coeff_Tspr = G_coeff_Tspr;
%         Migration.G_ldt = G_ldt;
%         Migration.exp = Params.exp;
%         Migration.time_const = time_const;
%         Migration.dec_rate = dec_rate;
%         Migration.inc_rate = inc_rate;
%         Migration.t_sim = Params.t_sim;
%         Migration.task_priorities = 1.8 + 0.7*rand(n_task,1);
%     end
% 
%     D = zeros(total_units,n_core);
%     P_d_tot = zeros(total_units,max_inst);
%     for c = 1:n_core
%         units_range = (c-1)*core_units+1:c*core_units;
%         D(units_range,c) = ones(core_units,1);
%         P_d_tot(units_range,:) = P_d(:,c,:);
%     end
%     
% %% Build parameter structure
% 
%     Scaling = struct(...
%         'inst', Params.inst);
%         
%     Power = struct(...
%         'P_d', P_d, ...
%         'P_d_tot', P_d_tot, ...
%         'P_s', P_s, ...
%         'vt', vt, ...
%         'K', K, ...
%         'D', sparse(D), ...
%         'T_range', Power.Temperature_range, ...
%         'vdd_tab', vdd_tab, ...
%         'Vdd_range', vdd_range', ...
%         'spd_range', spd_range', ...
%         'comp_method', 'core');        
%     
%     Sim_params = struct(...
%         'inst', Params.inst, ...
%         'IPC', Constant.IPC_SPEC, ...
%         'inst_sep', 10e6, ...
%         'n_core', n_core, ...
%         'n_task', n_task, ...
%         'core_units', core_units, ...
%         'total_units', total_units, ...
%         'policy', Params.policy, ...
%         'f_max', Params.f_max, ...
%         'vdd_max', vdd_max, ...
%         'T_amb', Thermal.t_amb_abs, ...
%         'T_ini', Params.T_ini, ...
%         'T_max', Thermal.t_max, ...
%         'spd_tol', 1/(10*Params.n_states), ...
%         'dtm_op', dtm_op, ...
%         't_step', t_step, ...
%         'method', Params.method, ...
%         'sub_method', sub_method, ...
%         'discrete', Params.discrete, ...
%         'n_spd_states', Params.n_states, ...
%         'ignore_leakage', false, ...
%         'spd_const', 1.8807e+03, ...
%         'A', Thermal.Hotspot.A, ...
%         'B', Thermal.Hotspot.B, ...        
%         'M', Thermal.Precomputed.M, ...
%         'R', Thermal.Precomputed.R, ...
%         'tau', Thermal.Time_constant.pkg, ...
%         'migration', Migration, ...
%         'scaling', Scaling, ...
%         'solver_options', optimset('Display','off'));
%     
%     % Special Cases
% 
%     if strcmp(Params.method, 'scaling') && strcmp(sub_method, 'accurate')
%         Power.P_d = zeros(total_units,1);
% 
%         for c = 1:n_core
%             units_range = (c-1)*core_units+1:c*core_units;
%             Power.P_d(units_range) = P_d(:,c);
%         end
%     elseif strcmp(sub_method, 'const_throttle')
%         A_HS_prime = (Sim_params.M - eye(total_units)) ...
%             * Thermal.Hotspot.B/Sim_params.R;
%         Sim_params.R_ss = (A_HS_prime\ (expm(A_HS_prime * 150) ...
%             - eye(total_units))) * Thermal.Hotspot.B;
%     end
%     
%     if strcmp(dtm_op, 'migration')
%         Power.comp_method = 'mig-total';
%     end
%     
%     function [err] = error_spd(vdd_in)
%         err = curr_spd - 1.6755e+03 * (vdd_in-vt)^1.2 / ...
%             (vdd_in*curr_temp^1.19);
%     end
% 


%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================






T_ini =Params.T_ini

P_dyn = Power1.P_d
 
%I_split =    %???
 
%C =    % computed by extact_hotspot_model
M = Precomputed.M
%S = 
%Thread_Location =
%Die =    % computed above
%Power =   % computed above
%Thermal =  % computed above
%Simulation =  % computed above
% 


% 
% T_ini = temperature_step_response_PWL_LDT(T_ini, P_dyn, I_split, ...
%              C, M, S, Thread_Location, Die, Power, Thermal, Simulation)
%          



% TEMPERATURE_STEP_RESPONSE_PWL_LDT computes the step response for the
% Hotspot thermal equivalent circuit over a specified time duration, for a
% given dynamic power vector and leakage-temperature piecewise-linear
% coefficients.
%
%   T_end = temperature_step_response_PWL_LDT(T_ini, P_dyn, I_split, ...
%             C, M, S, Thread_Location, Die, Power, Thermal, Simulation)
%   returns the temperature at the end of the time step whose value is
%   equal to Simulation.n_step_temp * Simulation.t_step_temp. The
%   temperature is actually computed at the intermediate n_step_temp steps
%   to reduce the error due to the piece-wise linear approximation of
%   leakage dependence on temperature (PWL-LDT). More on this below. 
%   The inputs to this function are:
%       T_ini - a column vector of size Thermal.N_units.total that
%       specifies the temperature of each Hotspot block at the beginning of
%       the time step.
%       P_dyn - a column vector of size Thermal.N_units.total that
%       specifies the dynamic power consumed by each Hotspot block when the
%       corresponding core operates at maximum speed. For non-die blocks,
%       the corresponding P_dyn elements are 0.
%       I_split, C, and M - column vectors of size Thermal.N_units.die that
%       together specify the PWL LDT regions and coefficients in each
%       region. For each die block i, I_split(i) denotes the PWL region 
%       corresponding to it initial temperature T_ini(i). The temperature
%       of the left boundary of this region is computed as T_left(i) =
%       Power.T_PWL_LDT(I_split(i)). The leakage power consumed by this 
%       block is given by 
%           P_leak(i) = C(i) + M(i) * (T_ini(i) - T_left(i))
%       S = a column vector of length Die.n_core that denotes the
%       (constant) speed of each core over the given time interval.
%       Thread_Location - a matrix of size Die.n_core_y x Die.n_core_x,
%       with Thread_Location(i, j) denoting the thread_id of the thread
%       running on core (i, j) over the given time interval. If the core is
%       inactive, Thread_Location(i, j) = 0.
%       Die, Power, Thermal, Simuation - Model parameter structures
%       generated by the function compute_model_parameters.m.
%   The temperature computation proceeds as follows. First, the given time
%   interval is divided into Simulation.n_step_temp steps, each of duration
%   Simulation.t_step_temp. The power consumption over time step j is
%   assumed to be
%       P_j(t) = diag(P_dyn) * D * S + C_{j-1} + diag(M_{j-1}) * (T_j(t) -
%       T_{left,j-1}), for all 0 <= t <= Simulation.t_step_temp
%   i.e. the PWL LDT coefficients for time step j are computed based on the
%   temperature at the beginning of the time step. This is an approximation
%   as the region in which T_j(t) lies could differ from that of T_{j-1}
%   (which is equal to T_j(0)).
%   The power and temperature vectors for time step j are related by
%       dT_j(t)/dt = A_HS * T_j(t) + B_HS * P_j(t), for all 0 <= t <=
%           Simulation.t_step_temp
%   where [A_HS, B_HS] represents the Hotspot thermal equivalent circuit.
%   Substituting for P_j in the above equation, we get
%       dT_j(t)/dt = A * T_j(t) + B * P, for all 0 <= t <=
%           Simulation.t_step_temp
%   where A = A_HS + B_HS * diag(M_{j-1}), B = B_HS, and 
%       P = diag(P_dyn) * D * S + C_{j-1} - diag(M_{j-1}) * T_{left,(j-1)}.
%   The solution for this equation is given by
%       T_j(t) = expm(A * t) * T_{j-1} + inv(A) * (expm(A * t) - I) * B * P
%           for all 0 <= t <= Simulation.t_step_temp
%   The larger the number of time steps, the less the distance between T_j
%   and T_{j-1}, and hence the more accurate our temperature estimation
%   that computes leakage power solely on T_{j-1}. But this comes at the
%   cost of greater simulation time.
%
% Ravishankar Rao, Arizona State University
% Created, Jan 31 2008

% function T_ini = temperature_step_response_PWL_LDT(T_ini, P_dyn, I_split, ...
%              C, M, S, Thread_Location, Die, Power, Thermal, Simulation)

% For each time step
for i_step = 1:Simulation.n_step %Simulation.n_step_temp
    % Compute matrices needed for temperature computation
    A = Thermal.Hotspot.A + Thermal.Hotspot.B * diag(M);
    T_left = Power.T_PWL_LDT(I_split);
    P = diag(P_dyn) * Thermal.Precomputed.D * S + C - diag(M) * T_left;
    expm_A_t = expm(A * Simulation.t_step_temp);
    % Find temperature at the end of the time step
    T_ini = expm_A_t * T_ini + ...
        inv(A) * (expm_A_t - Thermal.I) * Thermal.Hotspot.B * P;
    % Replicate this matrix
    T_ini_rep = repmat(T_ini, 1, Power.n_PWL_LDT+1);
    % Compare the above matrix with the PWL temperature regions to
    % determine the regions for the newly computed temperature vector
    I_split_new = sum(Power.T_PWL_LDT_rep <= T_ini_rep, 2);
    
    Block_index = Thermal.Block_index_core_1;
    for core_y = 1:Die.n_core_y
        for core_x = 1:Die.n_core_x
            Block_region = I_split_new(Block_index);
            % If any blocks crossed their previous PWL temperature regions, 
            % compute their new PWL LDT coefficients
            if any(Block_region ~= I_split(Block_index))
                thread_id = Thread_Location(core_y, core_x);  
                if thread_id
                    Thread_id = thread_id * Thread.Ones; 
                    Linear_index = ...
                        sub2ind(size(Thread.Power.ActLeak.C), ...
                            Block_region, Thermal.Block_index_core_1_t, ...
                                Thread_id);
                    C(Block_index) = Thread.Power.ActLeak.C(Linear_index);
                    M(Block_index) = Thread.Power.ActLeak.M(Linear_index);
                else
                    Linear_index = sub2ind(size(Power.Mean.StdLeak.C), ...
                        Block_region, Thermal.Block_index_core_1_t);
                    C(Block_index) = Power.Mean.StdLeak.C(Linear_index);
                    M(Block_index) = Power.Mean.StdLeak.M(Linear_index);
                end
                Block_index = Block_index + Thermal.N_units.per_core;
            end
        end
    end    
    I_split = I_split_new;
end
        
        
        
        


%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================

%--------------------------------------------------------------------------
% SCALING_EXP computes the optimal speeds and voltages of all cores to
% minimize the overall makespan of all tasks or equivalently maximize the
% overall throughput.
%
%   [Result] = scaling_exp(Sim_params, Power)
%
%   reads simulation and model parameters from 'Sim_params', and power
%   profile from 'Power' structures respectively. To learn more about the 
%   structures, type'help collect_sim_params.m'.
%
%   'Result' contains the speed, voltage, power and temperature values for
%   the simulation time.
%
% There are two options to compute the speeds, 'approximate' or 'accurate' 
% method. 
% 
% Read the comments in the file collect_sim_params.m.m to get a list
% and description of all model parameters computed. Changes to model
% parameters normally will not have to be made in
% collect_sim_params.m, but only in read_model_parameters.m. A
% typical usage would be:
%   Constant = read_model_parameters();
%   Constant.some_parameter = new_value;
%   [Die, Power, Thermal, Thread, Simulation] = ...
%        compute_model_parameters(Constant);
% 
% See start.m for a usage example.
%
% Vinay Hanumaiah, Arizona State University
%
% 
% 
% function [Result] = scaling_exp(Sim_params, Power)

%% Preliminaries 

    Method = Sim_params.method;
    Sub_method = Sim_params.sub_method;
    core_units = Sim_params.core_units;

    % Check method
    if strcmp(Method, 'dfs')
        method_str = 'DFS';
    elseif strcmp(Method, 'dvfs')
        method_str = 'DVFS';
    else
        error('Unidentified scaling specification in Record.sub_method');
    end

    % Check sub-method
    if strcmp(Sub_method, 'approx')
        sub_method_str = 'approximate';
        
        n_core = Sim_params.n_core;
        R = zeros(core_units,core_units,n_core);
        
        for c = 1:n_core
            die_units_i = (c-1)*core_units+1:c*core_units;
            R(:,:,c) = Sim_params.R(die_units_i,die_units_i);
        end
        Sim_params.scaling.R = R;
    elseif strcmp(Sub_method, 'accurate')
        sub_method_str = 'accurate';
    elseif strcmp(Sub_method, 'const_throttle')
        sub_method_str = 'const_throttle';
        Result = compute_ss_temp(Sim_params, Power);
        Sim_params.scaling.spd_ss = Result.spd;
    else
        error('Unidentified accuracy specification in Record.sub_method');
    end
    
    fprintf('Performing %s optimization using %s method...\n', method_str, ...
        sub_method_str);
    
    n_core = Sim_params.n_core;
    inst = Sim_params.scaling.inst;
    t_step = Sim_params.t_step;
    f_max = Sim_params.f_max;
    total_units = Sim_params.total_units;

    max_inst = max(inst);
    max_steps = 15*round(max(inst) / (f_max * t_step));
    
    P_d_tot = Power.P_d_tot;
    P_d = Power.P_d;
    Power.P_d_tot = zeros(total_units,1);
    Power.P_d = zeros(core_units,n_core);
    
    % Pre-assigning sizes
    spd = zeros(max_steps, n_core)+1;
    T = zeros(max_steps,total_units)+Sim_params.T_ini;
    P = zeros(max_steps,total_units);
    vdd = zeros(max_steps, n_core);
    
%% Simulation
    t = 0;
    dec_comp = 1;
    
    while any(inst > 0) 
        t = t+1;
        T_prev = T(t,:)';
        Sim_params.scaling.s_prev = spd(t,:)';
        
        if t > 1
            T_prev = T(t-1,:)';          
            Sim_params.scaling.s_prev = spd(t-1,:)';
        end
        
        inst_int = floor((Sim_params.inst - inst)./ ...
                (Sim_params.IPC * Sim_params.inst_sep)) + 1;
            
        T_diff = Sim_params.T_max - Sim_params.M*T_prev;
        
        for c = 1:n_core
            units_range = (c-1)*core_units+1:c*core_units;
            if inst_int(c) > size(P_d_tot,2)
                curr_inst_int = size(P_d_tot,2);
            else
                curr_inst_int = inst_int(c);
            end
            Power.P_d_tot(units_range) = P_d_tot(units_range,curr_inst_int);
        end
        
%         Sim_params.P_prev = P(t_prev,:)';
%         Sim_params.vdd_prev = vdd(t_prev,:);
%         Sim_params.spd_prev = spd(t_prev,:);

        if strcmp(Sim_params.sub_method, 'approx')
%             tic;
            for c = 1:n_core
                Power.P_d(:,c) = P_d(:,c,inst_int(c));
            end
            [spd(t,:),vdd(t,:),P(t,:),T(t,:)] = ...
                approx_dvfs(T_prev, T_diff, Sim_params, Power);
%             t_app = toc;
        elseif strcmp(Sim_params.sub_method, 'accurate')
%             tic
%             Sim_params.Power.comp_method = 'total';
%             Sim_params.scaling.T = 35*ones(total_units,1);
%             for i = 1:100
%                 spd = [i/100;i/100;i/100;i/100];
%                 [T,P,vdd] = get_power_temp(spd, Sim_params);
%                 obj(i) = sum(spd)*10-sum(P);
%             end
%             plot(obj);
%             return;
            
%% Some plots           
%             Sim_params.Power.comp_method = 'total';
% 
%             range = linspace(0.01,1,20);
%             range_len = length(range);
%             f = zeros(range_len,range_len);
%             
%             Sim_params.scaling.T = 45*ones(Sim_params.total_units);

%             range = linspace(0.01,1,100);
%             P_1 = zeros(range_len,1);
            
%             ind1 = 0;
%             for s1 = range
%                 ind1 = ind1 + 1
%                 spd = [s1;0];
%                 [T,P] = get_power_temp(spd, Sim_params);
%                 P_1(ind1) = sum(P);
%             end
% 
%             figure; plot(range,P_1)       
% %             
%             ind1 = 0;
%             for s1 = range
%                 ind1 = ind1 + 1
%                 ind2 = 0;
%                 for s2 = range
%                     ind2 = ind2 + 1;
%                     spd = [s1;s2];
%                     [T,P] = get_power_temp(spd, Sim_params);
%                     f(ind1,ind2) = sum(spd)/sum(P);
%                 end
%             end
% 
%             figure; surf(range,range,f)
            
%             f = zeros(range_len,1);
%             for T_in = [5 25 45 65]
%                 Sim_params.scaling.T = T_in*ones(Sim_params.total_units);
%                 ind = 0; T_in
%                 for spd = range
%                     ind = ind + 1;
%                     [T,P] = get_power_temp(spd, Sim_params);
%                     f(ind) = spd/sum(P);
%                 end
%                 plot(range, f); hold on;
%             end
%             legend('40','60','80','100');

%% back again            
%             tic
            [spd(t,:),vdd(t,:),P(t,:),T(t,:)] = ...
                accurate_dvfs(T_prev, Sim_params, Power);
%             t_acc = toc;
        end
%         fprintf('Improvement = %f\n',(t_acc-t_app)/t_app*100);
%         [max_diff, max_ind] = max(abs(spd_acc-spd_app));
%         spd_acc - spd_app
%         fprintf('Accuracy = %.5f\n', max_diff/spd_acc(max_ind)*100);        
%         
% return

        inst = inst - spd(t,:)'.*Sim_params.IPC*f_max*t_step;
        if any(spd(t,:) < -f_max*t_step)
            disp('check');
        end
        Sim_params.scaling.inst = inst;
        
        inst_comp = floor((max_inst - max(inst))/max_inst*10);
        
        if inst_comp >= dec_comp
            fprintf('%d%%\t', dec_comp*10);
            dec_comp = dec_comp +1;
        end
    end
    fprintf('\n');

    d_range = 1:t;
    Result = struct(...
        'spd', spd(d_range, :), ...
        'vdd', vdd(d_range, :), ...
        'P', P(d_range, :), ...
        'T', T(d_range, :));

%% Plot
    spd_volt_temp_plot(Result, Sim_params);

%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================



