function view_ptrace (ptrace_file, figno, scalefactor)
%View power tarce files data 
%ptrace_file = is the ptrace file containing thebechmark power
% figno = is the figure number where the data is displayed
% sclaefactor = a factor which scales the power and displays the scaled
% power
if nargin ==1
    
    figno=1;
    sf=1.0;
elseif nargin==2
    sf=1.0;
elseif nargin==3
    sf=scalefactor;
else 
    error ('Incorrect no of arguments');
end

cd('/Users/santanusarma/Dropbox/MATLAB/magma_v2/power_data/ptrace') ;


%check if the ptrace file is among the benchmark database

for i=1:26
    if strcmp(ptrace_file, id2ptrace(i))
       benchmark_found=1;
       break;
    else
       benchmark_found=0;
    end
    
end

%Read the trace file 

filename=ptrace_file;

%find the name and extension of the file name
%to be used to create the mat db file
dotpoint=strfind(filename, '.');
name=filename(1:dotpoint-1);
ext=filename(dotpoint+1:end);

%import the power data for the benchmark
power= importdata(filename, '\t', 1);
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

power1=(power.data).*sf;
[px, py]=size(power1)
power2=power1.';
avg_power=sum(sum(power2))/px;
avgpower=ones(1,px).*avg_power;


if benchmark_found
    %Read the average power from the *.p file for the benchmarks
    fp_pa= fopen([name, '.p'],'r');
    pow_avg=textscan(fp_pa, '%s%f');
    fclose(fp_pa);
    
    %Read the steady state power of the given benchmark
    power_ss= importdata([name,'.steady_ptrace'], '\t', 1);
    
    
    figure(figno)
    set (figno, 'Color', [1.0, 1.0, 1.0])
    subplot(2,2,1)
    mesh(power2)
    zlabel('Power in Watts')
    ylabel('Unit #')
    xlabel('time')
    axis tight
    subplot(2,2,2)
    plot(sum(power2))
    ylabel('Power in Watts')
    xlabel('time')
    title('Total Power over time')
    grid; axis tight
    hold on;
    plot(avgpower, '--r', 'LineWidth', 2)
    hold off;
    subplot(2,2,3)
    bar(sum(power1)/px)
    ylabel('Power in Watts')
    xlabel('Unit #')
    title('Average Unit Power over Time')
    axis tight
    subplot(2,2,4)
    bar(power_ss.data)
    xlabel('Unit #')
    ylabel('Power in Watts')
    title('Steady State Unit Power')
    axis tight
    
else
    
    figure(figno)
    set (figno, 'Color', [1.0, 1.0, 1.0])
    subplot(3,1,1)
    mesh(power2)
    zlabel('Power in Watts')
    ylabel('Unit #')
    xlabel('time')
    axis tight
    subplot(3,1,2)
    plot(sum(power2))
    ylabel('Power in Watts')
    xlabel('time')
    title('Total Power over time')
    grid; axis tight
    hold on;
    plot(avgpower, '--r', 'LineWidth', 2)
    hold off;
    subplot(3,1,3)
    bar(sum(power1)/px)
    ylabel('Power in Watts')
    xlabel('Unit #')
    title('Average Unit Power over Time')
    axis tight
    
end

end