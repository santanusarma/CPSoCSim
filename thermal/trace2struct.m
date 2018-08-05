function [struct, data]= trace2struct(trace_filename,trace_type,scalefactor)
%convert a trace file to struct data structure 

%scalefactor sf scales the power
if nargin ==1
    sf=1.0;
    trace_type='ttrace';
elseif nargin==2
    sf=1.0;
    if ~(strcmp(trace_type,'ptrace') || strcmp(trace_type,'ttrace'))
        error('invalid trace type')
    end

elseif nargin==3
    sf=scalefactor;
    if ~(strcmp(trace_type,'ptrace') || strcmp(trace_type,'ttrace'))
        error('invalid trace type')
    end
        
else
     error('invalid arguments')
end

if isa(trace_filename, 'char')
    
    %Read the trace file
    
    filename=trace_filename;
elseif isa(trace_filename, 'numeric')
    
    % power_traces is the id of the file
    if strcmp(trace_type, 'ptrace')
    filename= id2ptrace(trace_filename)
    elseif strcmp(trace_type, 'ttrace')
    filename= id2ttrace(trace_filename)
    else
       error('Invalid trace file id') 
    end

else
    error('Invalid trace filename')
end

    
    %find the name and extension of the file name
    %to be used to create the mat db file
    dotpoint=strfind(filename, '.');
    name=filename(1:dotpoint-1);
    ext=filename(dotpoint+1:end);
    
    %import the power data for the benchmark
    trace_struct= importdata(filename, '\t', 1)
    %trace_struct= importdata(filename, ' ', 1);
    if isfield(trace_struct, 'data')
        [px, py]=size(trace_struct.data);
        
    else
        trace_struct= importdata(filename, ' ', 1)
        if isfield(trace_struct, 'data')
            [px, py]=size(trace_struct.data);
            
        else
            error('Error in reading ptarce data');
        end
    end
    
    if strcmp(trace_type,'ptrace')
    %struct=trace_struct;
    data=(trace_struct.data)*sf;
    [px, py]=size(data);
    struct.power=data;
    struct.NoOfSamplingPoints=px;
    struct.NoOfBlocks=py;
    
    
    elseif strcmp(trace_type,'ttrace')
    %struct=trace_struct;
    data=(trace_struct.data)*sf;
    [px, py]=size(data);
    struct.temp=data;
    struct.NoOfSamplingPoints=px;
    struct.NoOfBlocks=py;
    
    
    else
        
     error('invalid trace type')   
    end
    
    if (py==length(trace_struct.textdata))
        struct.Labels=trace_struct.textdata;
    else
    struct.Labels=strsplit(trace_struct.textdata{:});
    end
    


end