%Thermal grid simulation of spec bechmarks for EV6 florplan 

%PowerTrace= benchmarks();


% %Power traces for the benchmarks
% PowerTrace={'ammp.ptrace'	'eon.ptrace'	'gcc.ptrace'	'parser.ptrace'	'vpr.ptrace'...
% 'applu.ptrace'	'equake.ptrace'	'gzip.ptrace'	'perlbmk.ptrace'	'wupwise.ptrace'...
% 'apsi.ptrace'	'facerec.ptrace'	'lucas.ptrace'	'sixtrack.ptrace'...
% 'art.ptrace'	'fma3d.ptrace'	'mcf.ptrace'	'swim.ptrace'...
% 'bzip2.ptrace'	'galgel.ptrace'	'mesa.ptrace'	'twolf.ptrace'...
% 'crafty.ptrace'	'gap.ptrace'	'mgrid.ptrace'	'vortex.ptrace', 'zeros.ptrace'}
% 
% for i=1:length(PowerTrace)
%     ptrace=['../power_data/ptrace/', PowerTrace{i}]
%     Simtime(i)=hotspot_grid_sim(ptrace)
% end
% 


%==========================================================================
% Thermal simulation of MPSOC 2x2 for all same benchmarks

nx=2
ny=2
NoOfBenchmarks=25
NoOfTraces =3000
[mpsoc_flp mpsoc_flp_file]=create_mpsoc_flp('ev6.flp',nx,ny)
mpsoc_config_file=select_config_file(mpsoc_flp_file,nx,ny)

for j=1: NoOfTraces
power_traces=randi(NoOfBenchmarks,nx)



%power_traces =i*ones(nx,ny)  %same benchmark for all

%mpsoc_config_file= 'hotspot_2x2_mpsoc.config'
%view_flp(mpsoc_flp_file)
% mpsoc_flp_filename='mpsoc_ev6_4x4.flp'
% cell2flp(mpsoc_flp, mpsoc_flp_filename)
%cd ('/Users/santanusarma/Dropbox/MATLAB/magma_v2/power_data/ptrace')
cd ('../power_data/ptrace')
[mpsoc_ptrace, powertr_file]=create_mpsoc_ptrace('ev6_hetero.flp',power_traces, nx,ny);
ss_trace=[getfilename(powertr_file),'.ss_ttrace']
%view_ptrace(powertr_file)

%cd ('/Users/santanusarma/Dropbox/MATLAB/magma_v2')
% temptr_file=thermalsim(mpsoc_flp_file, powertr_file, mpsoc_config_file);
% view_ttrace(temptr_file)

cd ('../../Hotspot')
ptrace=['../power_data/ptrace/', powertr_file]

%simtime=hotspot_grid_sim(ptrace,mpsoc_flp_file,mpsoc_config_file,getfilename(powertr_file),'64','64', '')
simtime=hotspot_grid_sim(ptrace,mpsoc_flp_file,mpsoc_config_file,ss_trace,'64','64', '')
history_pr.simtime{j}=simtime;
history_pr.powertr_file{j}=powertr_file;
history_pr.tempss_file{j}=ss_trace;
%view_grid_ttrace([getfilename(powertr_file),'.grid64x64.steady'], 64, 64, 'mesh', 'Centigrade')
eval(['!mv ',ss_trace, '  ../power_data/ptrace/mpsoc_traces/'])
cd ('../power_data/ptrace')
eval(['!mv ', getfilename(powertr_file),'.grid64x64.steady', '  ./mpsoc_traces/'])
eval(['!mv ', powertr_file, '  ./mpsoc_traces/'])
save('./mpsoc_traces/history_pr.mat',  'history_pr')
cd ('../../Hotspot')

end



