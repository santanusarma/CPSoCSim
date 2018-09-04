%Create grid thermal db

%Change to the thermal traces dir
cd ('../power_data/ptrace/mpsoc_traces/')
%Make a list of the grid trace file
!ls mpsoc_2x2*.grid32x32.steady > trace_list.txt
%read one-by-one each trace file


fd=fopen('trace_list.txt', 'r')

ThermalDb=[];
i=1;
while ~feof(fd)
    tline=fgetl(fd);
    if ~isempty(tline)
        grid_data=load(tline);
        Tss=grid_data(:,2);
        ThermalDb=[ThermalDb, Tss];
        ThermalGridDb.Tss{i}=Tss;
        ThermalGridDb.benchmarks{i}=tline;
        i=i+1;
    end
end
ThermalGridDb.ThermalDb=ThermalDb;
save('ThermalDb.mat', 'ThermalDb')
save('ThermalGridDb.mat', 'ThermalGridDb')
fclose(fd);
