%This is a test script for the function read_power_data.m


% The file format for the power data file. It consists of length(T)
% lines, with each line having the following format
%   DieTemperatureInC\tPDynAtMaxSpeed[Unit1]\tPActiveLeak[Unit1]
%   \tPStandbyLeak[Unit1]\tPDynAtMaxSpeed[Unit2]\tPActiveLeak[Unit2]
%   \tPStandbyLeak[Unit2] ... \tPDynAtMaxSpeed[UnitN]\tPActiveLeak[UnitN]
%   \tPStandbyLeak[UnitN]\n
% where N = n_units.


% Now lets read the power data file first from the file power.txt
filename='mcf.txt';
cd('power_data')
PowerData=load(filename)
cd ..
%PowerData = load('power.txt')

%Lets chech the size of the file
[x_rows y_cols]=size(PowerData)

% x_rows indicates the number of temperature poinst considered, while
% y_cols indicate the number of 3*blocks +1



% The first column of the data is the temperature points
TemperaturePoints= PowerData(:,1)

%The second col is the dynamic power consumed by the UNIT 1  and so on

% 1. If there are 20 units  y_coln should be 20*3 +1 =61
% 2. The dynamic power of each unit does not chnage with temperature, hence
% should be same under the second column.
% 3. Th epower file format for MiBench benchmarks need to be checked. The
% format for the Spec2000 only has been verified.

%Calculate the number of numits in the power file from the size and format of
%the file
n_units= (y_cols-1)/3

%Extract the Dynamic Power ( 2, 5, 8 .. elements) (P_d)
PowerDynamicIndx=  1 + 1:3:(3*n_units)
PowerDynamic=PowerData(1, PowerDynamicIndx)

%Extarct the Active Leakage Power /Static Power (P_a)
ActiveLeakagePowerIndx=  1 +PowerDynamicIndx
PowerActiveLeakage= PowerData(:, ActiveLeakagePowerIndx)

%Extarct the Standby Leakage Power /Inactive Power (P_s)
ActiveStandbyPowerIndx=  1 + ActiveLeakagePowerIndx
PowerStandbyLeakage= PowerData(:, ActiveStandbyPowerIndx)


%Get the names of the Units in the floor plan
 % Load the EV6 fllor plan 
[unit_name_1, width_1, height_1, x_left_1, y_bottom_1] = ...
    textread('ev6_1x1.flp', '%s %f %f %f %f');

Units=1:n_units;

%--------------------------------------------------------------------------
%generate the power trace files for hotspot simulation for the given floor
%plan
%--------------------------------------------------------------------------
%open a new file 
%find the name and extension of the file name
%hotspot_floor_plan='UltraSpacrT1.flp'
dotpoint=strfind(filename, '.');
name=filename(1:dotpoint-1);
ext=filename(dotpoint+1:end);
newfilename=[name, '.ptrace']
fp= fopen (newfilename, 'w')

%write the unit names in the first line
for i=1:n_units
tline=[unit_name_1{i}, ' '];
fwrite(fp,tline);
%fprintf(fp, tline)
end
fwrite(fp, char(10));
%write the power data for all the units line by line at each time stamp
dlmwrite(newfilename, PowerDynamic ,'-append', 'delimiter',' ') 

%--------------------------------------------------------------------------
%Generate trace file for the dynamic power for specified sampling points 
fp1= fopen ([name,'_dy','.ptrace'], 'w')
%write the unit names in the first line
for i=1:n_units
tline=[unit_name_1{i}, ' '];
fwrite(fp1,tline);
end
fwrite(fp1, char(10));
NoOfSamplingPoints=100;
for t=1:NoOfSamplingPoints
 for i=1:length(PowerDynamic)
    %dlmwrite(newfilename, PowerDynamic ,'-append', 'delimiter',' ') 
    tline=[num2str(PowerDynamic(i)), ' '];
    fwrite(fp1,tline);
 end
 fwrite(fp1, char(10));
end
%--------------------------------------------------------------------------







%--------------------------------------------------------------------------
%Plot the power data
%--------------------------------------------------------------------------
figure(1)
bar(PowerDynamic)
title('Dynamic Power')
xlabel('unit #')
ylabel('Power in Watts')
axis tight
set(gca, 'xtickLabel', unit_name_1) % check 



figure(2)
subplot(2,2,1)
bar(PowerDynamic)
title('Dynamic Power')
xlabel('Unit #')
ylabel('Power in Watts')
axis tight
%set(gca, 'xtickLabel', unit_name_1)

subplot(2,2,2)
surf(Units,TemperaturePoints, PowerActiveLeakage)
title('Active Lekage Power, Watts')
ylabel('Temperature T, in ^o C')
xlabel('Unit #')
zlabel('Active Lekage Power, Watts')
axis tight

subplot(2,2,3)
bar3(TemperaturePoints,PowerActiveLeakage)
title('Active Lekage Power, Watts')
ylabel('Temperature T, in ^o C')
xlabel('Unit #')
zlabel('Active Lekage Power, Watts')
axis tight


subplot(2,2,4)
surf(Units,TemperaturePoints, PowerStandbyLeakage)
title('Standby Lekage Power, Watts')
ylabel('Temperature T, in ^o C')
xlabel('Unit #')
zlabel('Standby Lekage Power, Watts')
axis tight

%--------------------------------------------------------------------------
%Display the power trace files as surface plots
%--------------------------------------------------------------------------

%Read the trace file 
cd('power_data/ptrace')
filename='gcc.ptrace'

%find the name and extension of the file name
%to be used to create the mat db file
dotpoint=strfind(filename, '.');
name=filename(1:dotpoint-1);
ext=filename(dotpoint+1:end);

%import the data 
power= importdata(filename, '\t', 1);

%Read the average power form the *.p file
fp_pa= fopen([name, '.p'],'r')
pow_avg=textscan(fp_pa, '%s%f')
fclose(fp_pa)

%Read the steady state power of the given benchmark
power_ss= importdata([name,'.steady_ptrace'], '\t', 1);

figure(3)
subplot(2,2,1)
mesh(power.data)
subplot(2,2,2)
surf(power.data)
subplot(2,2,3)
waterfall(power.data)
subplot(2,2,4)
ribbon(power.data)

power1=power.data;
[px, py]=size(power1)
power2=power1.';
avg_power=sum(sum(power2))/px;
avgpower=ones(1,px).*avg_power;

figure(4)
subplot(2,2,1)
mesh(power2)
zlabel('Power in Watts')
ylabel('Unit #')
xlabel('time')
axis tight
subplot(2,2,2)
plot(sum(power2))
ylabel('Power in Watts')
xlabel('time')
title('Total Power over time')
grid; axis tight
hold on;
plot(avgpower, '--r', 'LineWidth', 2)
hold off;
subplot(2,2,3)
bar(sum(power1)/px)
ylabel('Power in Watts')
xlabel('Unit #')
title('Average Unit Power over Time')
axis tight
subplot(2,2,4)
bar(power_ss.data)
xlabel('Unit #')
ylabel('Power in Watts')
title('Steady State Unit Power')
axis tight


% bar(pow_avg{2})
% xlabel('Unit #')
% ylabel('Power in Watts')
% title('Average Unit Power')
% axis tight


%--------------------------------------------------------------------------
%Thermal trace file gneration and ploting 

%Read the temperature trace file 
cd('/Users/santanusarma/Dropbox/MATLAB/magma_v2/power_data/ptrace')
filename='gcc.ttrace'

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
ylim([30, 90])
legend({'Maximum', 'Average', 'Minimum'}, 'Location', 'SouthEast')
hold off;
subplot(2,2,3)
bar(sum(tempr1)/tx)
ylabel('Temperature in ^o C')
xlabel('Unit #')
title('Average Unit Temperatures over Time')
axis tight
ylim([ 40 90])
subplot(2,2,4)
bar(tempr_ss.data(1:ty)-273.15)
xlabel('Unit #')
ylabel('Temperature in ^o C')
title('Steady State Unit Temperatures')
axis tight
ylim([ 40 90])

%--------------------------------------------------------------------------


%Read the temeprature trace file in  data structure / cell
cd('power_data/ptrace')
filename='gcc.ttrace'
%read the power trace file
fp_src=fopen(filename, 'r')
%Read the file into memory
%check for the first char # and empty line and ignore them
nheaderlines=1;
i=1;
while ~feof(fp_src) %ischar(tline)
     tline = fgets(fp_src);
     %disp(tline)
     
%if the line is not a comment then write to the file      
    
if ((tline(1)~='#') && (i>nheaderlines))  
    
    %The data in the file

    %disp(tline);
    %tline = fgets(fp_src);
    
    %Write the data without comments to the new file 
    numline=str2num(tline);
    data_cell{i}=numline;
    data_mat(i,:)=numline;
    
    %fwrite(fp, tline);
else
    %The comments line
    disp('Ignoring Header Lines & Comments')
    
end
i=i+1;

end
fclose(fp_src);
%--------------------------------------------------------------------------


