%Thermal grid simulation of spec bechmarks for EV6 florplan 

%PowerTrace= benchmarks();


%Power traces for the benchmarks
PowerTrace={'ammp.ptrace'	'eon.ptrace'	'gcc.ptrace'	'parser.ptrace'	'vpr.ptrace'...
'applu.ptrace'	'equake.ptrace'	'gzip.ptrace'	'perlbmk.ptrace'	'wupwise.ptrace'...
'apsi.ptrace'	'facerec.ptrace'	'lucas.ptrace'	'sixtrack.ptrace'...
'art.ptrace'	'fma3d.ptrace'	'mcf.ptrace'	'swim.ptrace'...
'bzip2.ptrace'	'galgel.ptrace'	'mesa.ptrace'	'twolf.ptrace'...
'crafty.ptrace'	'gap.ptrace'	'mgrid.ptrace'	'vortex.ptrace', 'zeros.ptrace'}

for i=1:length(PowerTrace)
    ptrace=['../power_data/ptrace/', PowerTrace{i}]
    Simtime(i)=hotspot_grid_sim(ptrace)
end