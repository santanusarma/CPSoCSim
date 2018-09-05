function [mpsoc_ptrace, powertr_file]=create_mpsoc_ptrace (single_core_flp, power_traces, nx,ny, ptrace_type)
% This function creates the mpsoc ptrace file and structure for given
% sequence of benchmarks and assignements to cores 

% Power_traces = a string contaning the file name of the benchmark e.g
% 'gcc.ptrace'
%              = a cell array containg a list of bencharks assigned to nx*ny
%             cores 
%              = a matrix contain the benchmark ids
% eg:     power_trace  ='gcc.ptrace'
%         power_traces ={ 'gcc.ptrace', 'gzip.ptrace';
%                         'zeros.ptrace', 'crafty.ptrace'}
%                      = {3 5; 0 23} % using ids 
%                      = [3 5; 0 23]
% single_core_flp = flp file/structure of the single core
% nx = no of cores in the x direction
% ny = no of cores in the y direction 

% ptrace_type= detremines if the trace is at 'block' level or 'core' level
% if 'core' level then teh power is summed up accross the block of the core


MinPowerPoints = 17414; %No of min point across all the traces

if nargin<4
    error('Incorrect Arguments')
elseif nargin==4
    ptracetype='block';
elseif nargin==5
    
    if strcmp(ptrace_type, 'block')
      ptracetype='block'  
    else
      ptracetype='core' 
    end
else
    error('Unknown Input Format')
end
    
%     power_traces='gcc.ptrace' % test
%     nx=2
%     ny=2
%     unit_name_N=[];

    [~, single_core_flp]= flp2struct(single_core_flp);
    if isa(power_traces, 'char')
        
        dotpoint=strfind(power_traces, '.');
        name=power_traces(1:dotpoint-1);
        ext=power_traces(dotpoint+1:end);
        
        path='../power_data/ptrace/';
        
        powertr_file= ['mpsoc_', num2str(nx),'x',num2str(ny),...
                        '_', name,'_all','.ptrace'];
        
        %import the power data for the benchmark
        power_tr= importdata(power_traces, '\t', 1);
        
        %repeat the power for nx*ny cores
        if strcmp(ptracetype, 'block')
            mpsoc_power=repmat(power_tr.data, 1, nx*ny);
        elseif strcmp(ptracetype,'core')
            mpsoc_power=repmat(sum(power_tr.data, 2), 1, nx*ny);
        else
            
            error('Unknwon ptracetype, other than block/core')
        end

        [sx, sy]=size(mpsoc_power);
        mpsoc_ptrace.power=mpsoc_power;
        
        %Generate the new lables
        if isa(single_core_flp, 'struct')
           
        elseif isa(single_core_flp, 'char') %flp file
            [~, single_core_flp]=flp2struct(single_core_flp)
            
        else
            error('Unknown Floorplan Input Format')
        end
        
        l=1;
        for i=1:nx
            for j=1:ny
                if strcmp(ptracetype, 'block')
                for k=1:single_core_flp.NoOfBlocks
                    unit_name_N{l} = [char(single_core_flp.Labels(k)),...
                        '_', num2str(i), '_', num2str(j)];
                    l=l+1;
                    %EmptyLabels{l}='';
                end
                elseif strcmp(ptracetype, 'core')
                    unit_name_N{l}=['core_', num2str(i), '_', num2str(j)];
                    l=l+1;
                else
                    error('Unknwon ptracetype, other than block/core')
                end
            end
        end
        
        mpsoc_ptrace.Labels = unit_name_N
        mpsoc_ptrace.benchmark=power_traces;
        
        mpsoc_ptrace.NoOfSamplingPoints=sx;
        mpsoc_ptrace.NoOfBlocks=sy;
        mpsoc_ptrace.nx=nx;
        mpsoc_ptrace.ny=ny;
        mpsoc_ptrace.filename=powertr_file;
        %write the data to file 
        pstruct2ptrace(mpsoc_ptrace, powertr_file);
        
    elseif isa(power_traces, 'cell') %cell containg the task assignment
        
        %In this format the power traces are assigned to the cores (
        %nx*ny)/ homogenous architecture using a cell array 
        
        % From the cell matrix the respective power trace is raed and the
        % pstruct data stucture is formed  and then the total power trace
        % file is generated 
        
        %Test conitions
        
%         power_traces ={ 'gcc.ptrace', 'gzip.ptrace';
%                         'zeros.ptrace', 'crafty.ptrace'}
%                      = {1 2; 5 8} % using ids 
%                     nx=2
%                     ny=2
%         
        [px, py]=size(power_traces);
        if (px~=nx) && (py~=ny)
            error('Power Trace Dimension Does not match the multicore floorplan')
        else
            
            l=1;
            power_db=[];
            ptrace_ids=[];
            for i=1:nx
                for j=1:ny
                    % Read the allocation of the traces
                    if isa(power_traces{i,j}, 'numeric')
                        powertr_filename=id2ptrace(power_traces{i,j});
                    else
                        powertr_filename=power_traces{i,j};
                    end
                    %ptrace_ids= strcat(ptrace_ids, num2str(ptrace2id(powertr_filename)))
                    
                    ptrace_ids=[ptrace_ids,'_', num2str(ptrace2id(powertr_filename))];
                    

                    %import the power data for the benchmark
                    power_tr= importdata(powertr_filename, '\t', 1);
                    if isfield(power_tr, 'data')
                        [px, py]=size(power_tr.data);
                    else
                        power_tr= importdata(powertr_filename, ' ', 1);
                        if isfield(power_tr, 'data')
                            [px, py]=size(power_tr.data);
                            
                        else
                            error('Error in reading ptarce data');
                        end
                    end
                    
                    if strcmp(ptracetype, 'block')
                        %Assign the power traces to the blocks of each core
                        x=power_tr.data;
                        power_db=[power_db, x(1:MinPowerPoints,:)];
                        
                        
                        for k=1:single_core_flp.NoOfBlocks
                            unit_name_N{l} = [char(single_core_flp.Labels(k)),...
                                '_', num2str(i), '_', num2str(j)];
                            l=l+1;
                            %EmptyLabels{l}='';
                        end
                        
                        
                    elseif strcmp(ptracetype, 'core')
                        x=sum(power_tr.data,2);
                        power_db=[power_db, x(1:MinPowerPoints,:)];
                        unit_name_N{l}=['core_', num2str(i), '_', num2str(j)];
                        l=l+1;
                        
                    else
                        error('Unknwon ptracetype, other than block/core')
                    end

                    
                end
            end
            
            [sx, sy]=size(power_db);
            
            mpsoc_ptrace.power=power_db;
            mpsoc_ptrace.Labels = unit_name_N;
            mpsoc_ptrace.benchmark=power_traces;
            
            mpsoc_ptrace.NoOfSamplingPoints=sx;
            mpsoc_ptrace.NoOfBlocks=sy;
            mpsoc_ptrace.nx=nx;
            mpsoc_ptrace.ny=ny;
            
            %The file name on which the dat has to be written
            if nx*ny<=8
                powertr_file=['mpsoc_',num2str(nx),'x',num2str(ny), ptrace_ids, '.ptrace']
            else
                powertr_file=['mpsoc_',num2str(nx),'x',num2str(ny),'_', num2str(sum(sum(power_traces))), '.ptrace']
            end
            pstruct2ptrace(mpsoc_ptrace, powertr_file);

        end
        
    elseif  isa(power_traces, 'numeric')  % Assignment made using the trace ids
        % eg power_trace =[1 2; 5 8] or 
        
        [px, py]=size(power_traces);
        if (px~=nx) && (py~=ny)
            error('Power Trace Dimension Does not match the multicore floorplan')
        else
            
            l=1;
            power_db=[];
            ptrace_ids=[];
            for i=1:nx
                for j=1:ny
                    % Read the allocation of the traces
                    powertr_filename=id2ptrace(power_traces(i,j))
                    
                    %ptrace_ids= strcat(ptrace_ids, num2str(ptrace2id(powertr_filename)))
                    
                    ptrace_ids=[ptrace_ids,'_', num2str(ptrace2id(powertr_filename))];
                    
                    
                    %import the power data for the benchmark
                    power_tr= importdata(powertr_filename, '\t', 1);
                    if isfield(power_tr, 'data')
                        [px, py]=size(power_tr.data);
                    else
                        power_tr= importdata(powertr_filename, ' ', 1);
                        if isfield(power_tr, 'data')
                            [px, py]=size(power_tr.data);
                            
                        else
                            error('Error in reading ptarce data');
                        end
                    end
                    
                    if strcmp(ptracetype, 'block')
                        x=power_tr.data;
                        power_db=[power_db, x(1:MinPowerPoints,:)];
                        
                        
                        
                        for k=1:single_core_flp.NoOfBlocks
                            unit_name_N{l} = [char(single_core_flp.Labels(k)),...
                                '_', num2str(i), '_', num2str(j)];
                            l=l+1;
                            %EmptyLabels{l}='';
                        end
                        
                    elseif strcmp(ptracetype, 'core')
                        x=sum(power_tr.data,2);
                        i, j
                        [sx, sy]=size(x)
                        power_db=[power_db, x(1:MinPowerPoints,:)];
                        unit_name_N{l}=['core_', num2str(i), '_', num2str(j)]
                        l=l+1;
                        
                    else
                        error('Unknwon ptracetype, other than block/core')
                    end
                    
                end
            end
            
            [sx, sy]=size(power_db);
            
            mpsoc_ptrace.power=power_db;
            mpsoc_ptrace.Labels = unit_name_N;
            mpsoc_ptrace.benchmark=power_traces;
            
            mpsoc_ptrace.NoOfSamplingPoints=sx;
            mpsoc_ptrace.NoOfBlocks=sy;
            mpsoc_ptrace.nx=nx;
            mpsoc_ptrace.ny=ny;
            
            %The file name on which the dat has to be written
            if nx*ny<=8
                powertr_file=['mpsoc_',num2str(nx),'x',num2str(ny), ptrace_ids, '.ptrace']
            else
                powertr_file=['mpsoc_',num2str(nx),'x',num2str(ny),'_', num2str(sum(sum(power_traces))), '.ptrace']
            end
            mpsoc_ptrace.filename=powertr_file;
            %write the data to file
            pstruct2ptrace(mpsoc_ptrace, powertr_file);
            
        end

        
        
    else
        error('Unknown Input Format')
    end



end