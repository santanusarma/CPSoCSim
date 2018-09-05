function [scaled_power, scaled_ptrace]= scalepower(ptrace_file, sf, out_ptrace_file)

% Sf = scalefactor 
% out_ptrace_file = nameof the output sclaed ptrace file 
% scaled_power = matrix representing the scaled power
% scaled_ptrace = filename where the sclaed power is stored
 %               = same as out_ptrace_file when explecity specified in
 %               input


%Read the trace file 

% filename=ptrace_file;
% 
% %find the name and extension of the file name
% %to be used to create the mat db file
% dotpoint=strfind(filename, '.');
% name=filename(1:dotpoint-1);
% ext=filename(dotpoint+1:end);
% 
% %import the power data for the benchmark
% power= importdata(filename, '\t', 1);
% if isfield(power, 'data')
%     [px, py]=size(power.data);
% else
%     power= importdata(filename, ' ', 1);
%     if isfield(power, 'data')
%         [px, py]=size(power.data);
%         
%     else
%         error('Error in reading ptarce data');
%     end
% end
% 
% power1=(power.data).*sf;
% [px, py]=size(power1)
% power2=power1.';
% avg_power=sum(sum(power2))/px;
% avgpower=ones(1,px).*avg_power;
% 
% %scaled_power=power1;
% scaled_power.data=power1;
% scaled_power.textdata=power.textdata;

if nargin==3
    scaled_ptrace=out_ptrace_file;
elseif nargin==2
    filename=ptrace_file;
    dotpoint=strfind(filename, '.');
    name=filename(1:dotpoint-1);
    ext=filename(dotpoint+1:end);

    scaled_ptrace=[name, '_', num2str(sf), '.', ext];
else
    error('invalid arguments');
end
    

%Write the vales to the file
pstruct=ptrace2pstruct(ptrace_file,sf);
scaled_power=pstruct;
pstruct2ptrace(scaled_power, scaled_ptrace);

end