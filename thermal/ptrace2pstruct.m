function [pstruct, power1]= ptrace2pstruct(power_traces,scalefactor)
%convert a ptrace file to pstruct data structure 

%scalefactor sf scales the power
if nargin ==1
    sf=1.0;
elseif nargin==2
    sf=scalefactor
else
     error('invalid arguments')
end

if isa(power_traces, 'char')
    
    %Read the trace file
    
    filename=power_traces;
    
elseif  isa(power_traces, 'numeric')
    
    % power_traces is the id of the file
    filename= id2ptrace(power_traces)
    
else
    error('Invalid ptrace file format')
end
    %find the name and extension of the file name
    %to be used to create the mat db file
    dotpoint=strfind(filename, '.');
    name=filename(1:dotpoint-1);
    ext=filename(dotpoint+1:end);
    
    %import the power data for the benchmark
    power= importdata(filename, '\t', 1);
    %power= importdata(filename, ' ', 1);
    if isfield(power, 'data')
        [px, py]=size(power.data);
    else
        power= importdata(filename, ' ', 1);
        if isfield(power, 'data')
            [px, py]=size(power.data);
            
        else
            error('Error in reading ptarce data');
        end
    end
    
    pstruct=power;
    power1=(power.data)*sf;
    [px, py]=size(power1);
    pstruct.power=power1;
    pstruct.NoOfSamplingPoints=px;
    pstruct.NoOfBlocks=py;
    pstruct.Labels=strsplit(pstruct.textdata{:});
    
%     power2=power1.';
%     avg_power=sum(sum(power2))/px;
%     avgpower=ones(1,px).*avg_power;
    
    



end