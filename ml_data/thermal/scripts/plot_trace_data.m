function [TraceData] = plot_trace_data(trace_file)
%PLOT_POWER_TRACE plots the power values for each block in the EV6
%processor given the input trace files

% INPUTS:
% Power_trace_file is file containg the power traces for ach blocks in the
% Hotspot format.

% OUTPUT
% Plots the value of the power for each block as function of time


% Get the format of the file 
% The first line gives the name of the bloacs 
% form second line onwards the power value of each bloakc is listed at the
% simulation time steps


% Get the name of teh blocks


% Get the data 
%power_trace_file='gcc.ptrace'
fid=fopen(trace_file, 'r')
fid_new=fopen('power_trace.txt', 'wb')
block_name_fid=fopen('block_name.txt','wb');
line=fgetl(fid)
line1=[line, '   '];
fprintf(block_name_fid, '%s\r',line); 
fclose(block_name_fid)


%collect the unit names ina cell array , BlockName
%find the blankspace in the functional unit  names header
spaces=isspace(line1);
FUName='';
[Str NoOfUnits]=sscanf(line, '%s');
j=1;
k=1;
for i=1:length(line1)
    
    if spaces(i)==0
        FUName(j)=line1(i)
        j=j+1
    elseif ((spaces(i)==1) && (spaces(i-1)==0)) 
        j=1
        BlockName{k}=FUName
        FUName='';
        k=k+1;
   
        %Save FUName
        
    end
    
end


%Gnerate the format of string data to be read by the function sscanf or
%textread
str='%s'
data_format=''
for i=1:NoOfUnits
data_format=[data_format, str];
end
format_data=['''',data_format,''''];



%Write the numeric data without the header FU names
while (~isempty(line))
    line=fgetl(fid);
    if line==-1
        break;
    end
    fprintf(fid_new, '%s\r',line); 
end

% Collect the data into a matrix

trace=dlmread('trace.txt')
[time_steps NoOfUnits]=size(trace)




% Plot the power data
for i=1:NoOfUnits
    figure(i)
    plot(power_trace(:,i))
    title(['trace data for   ', BlockName{i}])
    grid on;
    pause (2)
    close all
end

TraceData={BlockName, trace};

end

