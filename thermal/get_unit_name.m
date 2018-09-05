function [Units_name, Labels]=get_unit_name(filename)
%get the nameof the units/subsystems in the floor plan or tracefile 
% filename = thermal trace file 

% filename= 'mpsoc_2x2_0_11_0_0.ttrace'
% figno=1

cd('/Users/santanusarma/Dropbox/MATLAB/magma_v2/power_data/ptrace');
if nargin ==1
    
    figno=1;
    
    if strcmp(filename, '/') % path is given
        
    else
        cd('/Users/santanusarma/Dropbox/MATLAB/magma_v2/power_data/ptrace');
    end
    
else
    
end

%cd('/Users/santanusarma/Dropbox/MATLAB/magma_v2/power_data/ptrace')
%filename='gcc_2x2.ttrace'; %'gcc.ttrace'

%find the name and extension of the file name
%to be used to create the mat db file
dotpoint=strfind(filename, '.');
name=filename(1:dotpoint-1);
ext=filename(dotpoint+1:end);

%import the temperature trace data 
tempr= importdata(filename, '\t', 1);

%Units_name = strsplit(tempr.)

%get the unit names from the trace file
Units_name = tempr.textdata;
Labels = LegendLabels(Units_name);

end