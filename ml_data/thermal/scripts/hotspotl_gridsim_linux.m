function [simtime] =hotspotl_gridsim_linux(ptrace_file, flp_file, config_file, steady_file, grid_rows, grid_cols, out_path)
% This function performs hotspot based grid simulation
% 
% COMMAND FORMAT
% !./hotspot -c hotspot.config -f ev6.flp -p gcc.ptrace -steady_file gcc.steady -model_type grid -grid_steady_file gcc.grid55x55.steady -grid_rows 55 -grid_cols 55
% !./grid_thermal_map.pl ev6.flp gcc.grid55x55.steady > gcc2_55x55.svg
% !./grid_thermal_map_v5p0.pl ev6.flp gcc.grid55x55.steady > gcc2_55x55.svg
% 

t_start = tic;

[filename, ext]= getfilename(ptrace_file)

% pathstr=strfind(ptrace_file, '\');
% dotpoint=strfind(ptrace_file, '.');
% if isempty(pathstr) && length(dotpoint)==1
%     filename=ptrace_file(1:dotpoint-1);
%     ext=ptrace_file(dotpoint+1:end);
% else
%     filename=ptrace_file(pathstr(end)+1:dotpoint(end)-1);
%     ext=ptrace_file(dotpoint(end)+1:end);
%     
% end

if nargin==1
    
    cmd_part1=[' -c hotspot.config ']
    cmd_part2=[' -f ev6.flp ']
    cmd_part3=[' -p ', ptrace_file]
    cmd_part4=[' -steady_file ', filename, '.steady ']
    cmd_part5=[' -model_type grid ']
    %cmd_part6=[' -grid_steady_file ', filename, '.grid64x64.steady ' ]
    file_w_path=['../power_data/ptrace/', filename]
    cmd_part6=[' -grid_steady_file ', file_w_path, '.grid64x64.steady ' ]
    cmd_part7=[' -grid_rows 64 -grid_cols 64 ']
    
elseif nargin==2
    cmd_part1=[' -c hotspot.config ']
    cmd_part2=[' -f ', flp_file]
    cmd_part3=[' -p ', ptrace_file]
    cmd_part4=[' -steady_file ', filename, '.steady ']
    cmd_part5=[' -model_type grid ']
    cmd_part6=[' -grid_steady_file ', filename, '.grid64x64.steady ' ]
    cmd_part7=[' -grid_rows 64 -grid_cols 64 ']
    
elseif nargin==3
    cmd_part1=[' -c ', config_file]
    cmd_part2=[' -f ', flp_file]
    cmd_part3=[' -p ', ptrace_file]
    cmd_part4=[' -steady_file ', filename, '.steady ']
    cmd_part5=[' -model_type grid ']
    cmd_part6=[' -grid_steady_file ', filename, '.grid64x64.steady ' ]
    cmd_part7=[' -grid_rows 64 -grid_cols 64 ']
elseif nargin==4
    cmd_part1=[' -c ', config_file]
    cmd_part2=[' -f ', flp_file]
    cmd_part3=[' -p ', ptrace_file]
    cmd_part4=[' -steady_file ', steady_file]
    cmd_part5=[' -model_type grid ']
    cmd_part6=[' -grid_steady_file ', filename, '.grid64x64.steady ' ]
    cmd_part7=[' -grid_rows 64 -grid_cols 64 ']
    
elseif nargin==5 %(check this case)
    cmd_part1=[' -c ', config_file]
    cmd_part2=[' -f ', flp_file]
    cmd_part3=[' -p ', ptrace_file]
    cmd_part4=[' -steady_file gcc.steady ']
    cmd_part5=[' -model_type grid ']
    
    grid_cols=steady_file;
    cmd_part6=[' -grid_steady_file ', filename, '.grid',grid_rows, 'x',grid_cols,'.steady ' ]
    cmd_part7=[' -grid_rows ', grid_rows,  ' -grid_cols ', grid_cols ]

    
elseif nargin==6
    cmd_part1=[' -c ', config_file]
    cmd_part2=[' -f ', flp_file]
    cmd_part3=[' -p ', ptrace_file]
    cmd_part4=[' -steady_file ', steady_file]
    cmd_part5=[' -model_type grid ']
    cmd_part6=[' -grid_steady_file ', filename, '.grid',grid_rows, 'x',grid_cols,'.steady ' ]
    cmd_part7=[' -grid_rows ', grid_rows,  ' -grid_cols ', grid_cols ]

elseif nargin==7
    cmd_part1=[' -c ', config_file]
    cmd_part2=[' -f ', flp_file]
    cmd_part3=[' -p ', ptrace_file]
    cmd_part4=[' -steady_file ', steady_file]
    cmd_part5=[' -model_type grid ']
    
    %filename='test.txt'
    %out_path='../power_data/ptrace/'
    if ~isempty(strfind(out_path, '/')) && exist(out_path,'dir')
        file_w_path=[out_path, filename]
        %disp('true path')
    else
        file_w_path=['../power_data/ptrace/', filename]
    end
    cmd_part6=[' -grid_steady_file ', file_w_path, '.grid', grid_rows, 'x',grid_cols,'.steady ' ]
    cmd_part7=[' -grid_rows ', grid_rows,  ' -grid_cols ', grid_cols ]

else
    error('Incorrect format')
    
end
cmd_part0='./hotspot_lin ';
cmd=[cmd_part0, cmd_part1, cmd_part2, cmd_part3, cmd_part4, cmd_part5, cmd_part6, cmd_part7]
system(cmd);


simtime=toc(t_start);

end