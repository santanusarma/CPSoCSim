function [temptr_fl]=thermalsim(floorplan, ptrace_file, config_file)
%Perform hotspot thermal simulation and produce the thermal simulation
%files 

% floorplan = floorplan file of the architecture
% ptrace_flie = ptrace file of the benchmark
% both the files should match their dimensions

if nargin==2
    
    config_file = 'hotspot.config'
elseif nargin==3
    
    if ~isa(config_file, 'char')
        error('Invalid Config File');
    end
    %config_file = 'hotspot.config'
else
    error('Invalid Number of Inputs');
end

if ~isa(floorplan, 'char')
    error('Invalid Floorplan File');
elseif ~isa(ptrace_file, 'char')
    error('Invalid Power Trace File');
end

%filename='mpsoc_2x2_3_8_0_23.ptrace'
filename=ptrace_file;
dotpoint=strfind(filename, '.');
name=filename(1:dotpoint-1);
ext=filename(dotpoint+1:end);

temptr_fl=[name, '.ttrace'];
powertr_file=['../power_data/ptrace/',filename];
temptr_file=['../power_data/ptrace/', name, '.ttrace'];
tempss_file=['../power_data/ptrace/', name, '.ss_ttrace'];
flp_file=['../power_data/ptrace/', floorplan];

cmd_part1= ['!./hotspot -c ', config_file, ' -f ../', floorplan];
cmd_part2= [' -p ', powertr_file];
cmd_part3= [' -o ', temptr_file];
cmd_part4= [' > ', tempss_file];

cmd=[cmd_part1, cmd_part2, cmd_part3, cmd_part4];

cd ('../hotspot/HotSpot-5.02');
disp('Performing Thermal Simulations ! Please wait for a moment')
tic;
t1=cputime;
eval(cmd);

%cd ('CPSOCSIM_ROOT')
t2=cputime;
simtime=t2-t1;
toc;
disp(['Thermal Simulations Completed in ', num2str(simtime), ' sec'])

end
