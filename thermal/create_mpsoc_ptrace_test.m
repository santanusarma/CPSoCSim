%Test of create_mpsoc_ptrace

[mpsoc_ptrace, powertr_file]=create_mpsoc_ptrace('ev6_hetero.flp','gcc.ptrace', 8, 8, 'core')

power_traces=randi(25, 8, 8)
%[px, py]=size(power_traces{:})

[mpsoc_ptrace, powertr_file]=create_mpsoc_ptrace('ev6_hetero.flp',power_traces, 8, 8, 'core')