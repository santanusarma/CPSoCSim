%HMP floorplan 

%Given the floorplan of two processor types and total area , construct a
%HMP floorplan.

%We have only one floor plan and the corresponding power traces. Th eother
%floor plan is created by scaling the EV6 floorplan by the approproate are
%and power scaling constants


%--------------------------------------------------------------------------
%Alpha Processor Relative Perfromance at 65 nm
%Alpha Core    Peak Power (W)   Average Power (W)   IPC     Area    Power
%EV4            4.97            3.73                1.00    1.00    1.00
%EV5            9.83            6.88                1.30    1.76    1.84
%EV6            17.8            10.68               1.87    8.54    2.86
%EV8            92.88           46.44               2.14    82.20   12.45

 AlphaDB= [ 4.97            3.73                1.00    1.00    1.00
            9.83            6.88                1.30    1.76    1.84
            17.8            10.68               1.87    8.54    2.86
            92.88           46.44               2.14    82.20   12.45];
%--------------------------------------------------------------------------
 
%Read the Ev6 and Ev4 floor plan with 18 blocks that have power traces 

cd ('power_data/ptrace')
ev6_hetero=flp2cell('ev6_hetero.flp')
%ev4_hetero=scalefloorplan('ev6_hetero.flp', 1/8.54)
ev4_hetero=scalefloorplan('ev6_hetero.flp', 1/3) % 3~sqrt(8.54)
cell2flp(ev4_hetero, 'ev4_hetero.flp')

[ev6_cell_flp ev6_str_flp]=flp2cell('ev6_hetero.flp')

view_flp('ev6_hetero.flp','',1)

view_flp('ev4_hetero.flp','',2)

cd ..
cd ..


%generate the homogenous multi-core floorplan of big core

create_floorplan ('power_data/ptrace/ev6_hetero.flpp', 2,2, 'mpsoc_ev6_2x2.flp')
create_floorplan ('power_data/ptrace/ev6_hetero.flpp', 2,1, 'mpsoc_ev6_2x1.flp')
create_floorplan ('power_data/ptrace/ev6_hetero.flpp', 1,2, 'mpsoc_ev6_1x2.flp')

%Gnerate the homogenous multi-core floorplan of the small core
create_floorplan ('power_data/ptrace/ev4_hetero.flpp', 3,3, 'mpsoc_ev4_3x3.flp')
create_floorplan ('power_data/ptrace/ev4_hetero.flpp', 3,6, 'mpsoc_ev4_3x6.flp')
create_floorplan ('power_data/ptrace/ev4_hetero.flpp', 6,3, 'mpsoc_ev4_6x3.flp')

%create the floorplan structures 
[ev6_cell, ev6_str]=flp2cell('power_data/ptrace/ev6_hetero.flp')
[ev4_cell, ev4_str]=flp2cell('power_data/ptrace/ev4_hetero.flp')
[ev4_cell1, ev4_str1]= flp_oper(ev4_str, 'shift', [2, 2], 'scale', 2)

ev6_2x2=create_mpsoc_flp(ev6_str, 2,2)
ev6_2x1=create_mpsoc_flp(ev6_str, 2,1)
ev6_1x2=create_mpsoc_flp(ev6_str, 1,2)
ev4_3x3=create_mpsoc_flp(ev4_str, 3,3)
ev4_3x6=create_mpsoc_flp(ev4_str, 3,6)
ev4_6x3=create_mpsoc_flp(ev4_str, 6,3)

view_flp_struct(ev4_3x3)


%Mix both the floor plans to form big-little
hmp_1x2_l=flpjoin(ev6_str_flp, ev4_3x3, 'left')
view_flp_struct(hmp_1x2_l)

hmp_1x2_r=flpjoin(ev6_str_flp, ev4_3x3, 'right')
view_flp_struct(hmp_1x2_r)

hmp_1x2_u=flpjoin(ev6_str_flp, ev4_3x3, 'up')
view_flp_struct(hmp_1x2_u)

hmp_1x2_d=flpjoin(ev6_str_flp, ev4_3x3, 'down')
view_flp_struct(hmp_1x2_d)


hmp_2x2_lr1=flpjoin(hmp_1x2_l, hmp_1x2_r, 'up')
view_flp_struct(hmp_2x2_lr1)

hmp_2x2_lr2=flpjoin(hmp_1x2_l, hmp_1x2_r, 'down')
view_flp_struct(hmp_2x2_lr2)

% used in cpsoc
hmp_2x2_ud1=flpjoin(hmp_1x2_u, hmp_1x2_d, 'right')
view_flp_struct(hmp_2x2_ud1)
cell2flp(hmp_2x2_ud1, 'hmp_2x2_ud1.flp')
cell2flp(hmp_2x2_ud1, 'cpsoc_2x2_ud1.flp')
view_flp('cpsoc_2x2_ud1.flp')

hmp_2x2_ud2=flpjoin(hmp_1x2_u, hmp_1x2_d, 'left')
view_flp_struct(hmp_2x2_ud2)

%generate the power trace for the multi-core
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%Perform thermal sim of homogenous big -core floorplan
filename='gcc_2x2.ptrace'
dotpoint=strfind(filename, '.');
name=filename(1:dotpoint-1);
ext=filename(dotpoint+1:end);

powertr_file= ['../power_data/ptrace/','gcc_2x2.ptrace']
temptr_file=['../power_data/ptrace/', name, '.ttrace']
tempss_file=['../power_data/ptrace/', name, '.ss_ttrace']
flp_file=['../power_data/ptrace/', 'mpsoc_ev6_2x2.flp']

cmd_part1= ['!./hotspot -c hotspot.config -f ../mpsoc_ev6_2x2.flp']
cmd_part2= [' -p ', powertr_file]
cmd_part3= [' -o ', temptr_file]
cmd_part4= [' > ', tempss_file]

cmd=[cmd_part1, cmd_part2, cmd_part3, cmd_part4]

cd ('/Users/santanusarma/Dropbox/MATLAB/magma_v2/HotSpot-5.02')
eval(cmd)

cd ('/Users/santanusarma/Dropbox/MATLAB/magma_v2')

%--------------------------------------------------------------------------
%plot the thermal traces 
%--------------------------------------------------------------------------
%Read the temperature trace file 
cd('/Users/santanusarma/Dropbox/MATLAB/magma_v2/power_data/ptrace')
filename='gcc_2x2.ttrace'; %'gcc.ttrace'

%find the name and extension of the file name
%to be used to create the mat db file
dotpoint=strfind(filename, '.');
name=filename(1:dotpoint-1);
ext=filename(dotpoint+1:end);

%import the temperature trace data 
tempr= importdata(filename, '\t', 1);

%Read the steady state power of the given benchmark
tempr_ss= importdata([name,'.ss_ttrace'], '\t', 1);

tempr1=tempr.data;
[tx, ty]=size(tempr1);
tempr2=tempr1.';

avg_tempr=sum(sum(tempr1)/tx)/ty;
avgtempr=ones(1,tx).*avg_tempr;

Tmin=min(min(tempr2))-5;
Tmax=max(max(tempr2))+5;

figure(5)
subplot(2,2,1)
mesh(tempr2)
zlabel('Temperature in ^o C')
ylabel('Unit #')
xlabel('time')
title(['Thermal Profile with ', name, ' benchmark'])
axis tight
subplot(2,2,2)
plot(max(tempr2))
ylabel('Temperature in ^o C')
xlabel('time')
title('Maximum Temperature over time')
grid; %axis tight
hold on;
% hold on;
plot(avgtempr, '--r', 'LineWidth', 2)
plot(min(tempr2), '--g', 'LineWidth', 2)
xlim([0 2.2e4])
ylim([30, Tmax])
legend({'Maximum', 'Average', 'Minimum'}, 'Location', 'SouthEast')
hold off;
subplot(2,2,3)
bar(sum(tempr1)/tx)
ylabel('Temperature in ^o C')
xlabel('Unit #')
title('Average Unit Temperatures over Time')
axis tight
ylim([ 40 Tmax])
subplot(2,2,4)
bar(tempr_ss.data(1:ty)-273.15)
xlabel('Unit #')
ylabel('Temperature in ^o C')
title('Steady State Unit Temperatures')
axis tight
ylim([ 40 Tmax])

%--------------------------------------------------------------------------
%Perform thermal sim of homogeneous small-core floorplan


%Perform theraml simulation for big-little architecture 


%Compute MTTF 



%Compute change in Vth and Aging 



%


