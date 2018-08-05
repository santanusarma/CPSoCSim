function  [A_HS, B_HS, G, Cap, Sys, SysD]=thermal_model(flp, config_file, ptrace_file)
%This function only works with hotspot v5.0.0

%INPUTS:
% flp = floor plan  file inhotspot format
% config_file= configuration file of hotspot

%OUTPUTS:
% A_HS = system coupling matrix 
% B_HS = the control matrix 
% G  = The conductance at each node of the bloack model
% Cap = capacitance at each node of the 
% Sys = state-space representaion of the thermal model in continuous time
% SysD = SS representaion in digital domain for given sampling time Ts

% Extract the A and B matrices for the state space representation of the
% RC thermal equivalent circuit. Hotspot uses the following state-space
% equation: dT/dt + inv(A) * B * T = inv(A) * P, 
% or dT/dt = -inv(A) * B * T + inv(A) * P. Comparing with the standard
% state space notation: dT/dt = A_HS * T + B_HS * P, we have
% A_HS = -inv(A) * B, B_HS = inv(A). The inv(A) and B matrices are computed
% by Hotspot and dumped into the files matrix_inva and matrix_b,
% respectively. 

%According to hotspot 5.0.2 convensions :
% /* creates 2 matrices: invA, C: dT + A^-1*BT = A^-1*Power, 
%  * C = A^-1 * B. note that A is a diagonal matrix (no lateral
%  * capacitances. all capacitances are to ground). also note that
%  * it is stored as a 1-d vector. so, for computing the inverse, 
%  * inva[i] = 1/a[i] is just enough. 
%  */

%From this A_HS =-inv(A)*B=C, and B_HS = inv(A) in hotspot convension

%==========================================================================

%INPUT PROCERSSING
if nargin==2 % the hotspot config file is provided
    cd ./HotSpot-5.02-mod
    if exist('model_vec_inva.txt')==2
    !mv model_vec_inva.txt model_vec_inva.txt.old
    !mv model_mat_b.txt model_mat_b.txt.old
    !mv model_vec_gamb.txt model_vec_gamb.txt.old
    !mv model_vec_a.txt model_vec_a.txt.old
    end
    cmd =['./thermal_model', ' -c ', config_file, ' -f ', flp];
    % The  input files dimension must match that of the floor plan
    if isa(flp, 'char')
    [~, flp_struct]=flp2struct(['../../', flp])
    
    %generate a dummp ptrace file of the same dimensions
    NoOfBlocks= flp_struct.NoOfBlocks;
    power_traces=randi(1, NoOfBlocks, 1)
    [mpsoc_ptrace, powertr_file]=create_mpsoc_ptrace(flp,power_traces, NoOfBlocks, 1, 'core')
    
    [name, etx]=getfilename(powertr_file)
    cmd =['./thermal_model', ' -c ', config_file, ' -f ', flp];
    % The  input files dimension must match that of the floor plan
   
    dummy_cmd=[' -p ', ptrace_file, ' -o ', name, '.ttrace', ' -steady_file ', name, '.steady'];
    cmd=[cmd,dummy_cmd];
    system(cmd);
    
    else %default case
    dummy_cmd=[' -p gcc.ptrace -o gcc.ttrace -steady_file gcc.steady'];
    cmd=[cmd,dummy_cmd];
    system(cmd);
    end
    
elseif nargin ==3
    cd ./HotSpot-5.02-mod
    [name, etx]=getfilename(ptrace_file)
    cmd =['./thermal_model', ' -c ', config_file, ' -f ', flp];
    % The  input files dimension must match that of the floor plan
   
    dummy_cmd=[' -p ', ptrace_file, ' -o ', name, '.ttrace', ' -steady_file ', name, '.steady']
    cmd=[cmd,dummy_cmd];
    system(cmd);
    
else
    error('Incorrect Arguments')
end




A_1_inv_vec = load('model_vec_inva.txt');
A_1_inv = diag(A_1_inv_vec);
B_1  = load('model_mat_b.txt');
A_HS = -A_1_inv*B_1;
B_HS = A_1_inv;

%Theconductance at each node of model is:
G= load('model_vec_gamb.txt');

%The capacitance at each node is stored in a[i] of the model
Cap= load('model_vec_a.txt');


[xa,ya]=size(A_HS) % n_units_in_die*2 +14
[xb yb]=size(B_HS) % n_units_in_die*2 +14

NoOfStates=xa;
% No of Extra Node = 12 in hotspot 5.0 and 4 layers are used to model the
% processor 

ExtraNodes=(NoOfStates-4*30)  

n_units_in_die=(NoOfStates -ExtraNodes)/4

C_HS=eye(NoOfStates);
D_HS=zeros(yb,NoOfStates);
x0=ones(NoOfStates,1)*(333.15-273.15);

%Write the system matrices to file
dlmwrite('A_HS.txt', A_HS)
dlmwrite('B_HS.txt', B_HS)
dlmwrite('C_HS.txt', C_HS)
dlmwrite('D_HS.txt', D_HS)
dlmwrite('X0.txt', x0)

%State-space representation of the thermal model (continuous time)
Sys =ss(A_HS, B_HS, C_HS, D_HS)

%Discretize the model;
Ts=3.333e-06  %sampling time of the discretization model, 
%to be read form config_file
%Ts=1.0e-4
SysD=c2d(Sys,Ts,'zoh') 
[A_HZ,B_HZ,C_HZ,D_HZ] = ssdata(SysD)
%Write the system matrices to file
dlmwrite('A_HZ.txt', A_HZ);
dlmwrite('B_HZ.txt', B_HZ);
dlmwrite('C_HZ.txt', C_HZ);
dlmwrite('D_HZ.txt', D_HZ);

cd ..
end