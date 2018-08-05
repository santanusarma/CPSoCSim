function ptrace_id= ptrace2id(filename)
%Given the file name, the function provides the id of teh benchmark

PowerTrace={'zeros.ptrace','ammp.ptrace'	'eon.ptrace'	'gcc.ptrace'	'parser.ptrace'	'vpr.ptrace'...
'applu.ptrace'	'equake.ptrace'	'gzip.ptrace'	'perlbmk.ptrace'	'wupwise.ptrace'...
'apsi.ptrace'	'facerec.ptrace'	'lucas.ptrace'	'sixtrack.ptrace'...
'art.ptrace'	'fma3d.ptrace'	'mcf.ptrace'	'swim.ptrace'...
'bzip2.ptrace'	'galgel.ptrace'	'mesa.ptrace'	'twolf.ptrace'...
'crafty.ptrace'	'gap.ptrace'	'mgrid.ptrace'	'vortex.ptrace'};

for i=1:length(PowerTrace)
    
    if strcmp(PowerTrace{i},filename)
        ptrace_id=i-1;
        break;
    else
        ptrace_id=0;
    end
    
end

end
