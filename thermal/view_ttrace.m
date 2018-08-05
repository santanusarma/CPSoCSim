function view_ttrace(filename, figno)
%View temperature tarce files data 

cd('/Users/santanusarma/Dropbox/MATLAB/magma_v2/power_data/ptrace');
if nargin ==1
    
    figno=1;
    
    if strcmp(filename, '/') % path is given
        
    else
        cd('/Users/santanusarma/Dropbox/MATLAB/magma_v2/power_data/ptrace');
    end
    
else
    
end

%cd('/Users/santanusarma/Dropbox/MATLAB/magma_v2/power_data/ptrace')
%filename='gcc_2x2.ttrace'; %'gcc.ttrace'

%find the name and extension of the file name
%to be used to create the mat db file
dotpoint=strfind(filename, '.');
name=filename(1:dotpoint-1);
ext=filename(dotpoint+1:end);

%import the temperature trace data 
tempr= importdata(filename, '\t', 1);

%Read the steady state power of the given benchmark
tempr_ss= importdata([name,'.ss_ttrace'], '\t', 1);

tempr1=tempr.data;
[tx, ty]=size(tempr1);
tempr2=tempr1.';

avg_tempr=sum(sum(tempr1)/tx)/ty;
avgtempr=ones(1,tx).*avg_tempr;

Tmin=min(min(tempr2))-5;
Tmax=max(max(tempr2))+5;


figure(figno)
set (figno, 'Color', [1.0, 1.0, 1.0])
subplot(2,2,1)
mesh(tempr2)
zlabel('Temperature in ^o C')
ylabel('Unit #')
xlabel('time')
title(['Thermal Profile with ', name, ' benchmark'])
axis tight
subplot(2,2,2)
plot(max(tempr2))
ylabel('Temperature in ^o C')
xlabel('time')
title('Maximum Temperature over time')
grid; %axis tight
hold on;
% hold on;
plot(avgtempr, '--r', 'LineWidth', 2)
plot(min(tempr2), '--g', 'LineWidth', 2)
xlim([0 tx])
ylim([30, Tmax])
legend({'Maximum', 'Average', 'Minimum'}, 'Location', 'SouthEast')
hold off;
subplot(2,2,3)
Tempr=sum(tempr1)/tx;
bar(Tempr)
ylabel('Temperature in ^o C')
xlabel('Unit #')
title('Average Unit Temperatures over Time')
axis tight
ylim([ 40 Tmax])
subplot(2,2,4)
Tempr_ss=tempr_ss.data(1:ty)-273.15;
bar(Tempr_ss)
xlabel('Unit #')
ylabel('Temperature in ^o C')
title('Steady State Unit Temperatures')
axis tight
ylim([ 40 Tmax])


%Plot the reliability / MTTF of the blocks
figure(figno+1)
set (figno+1, 'Color', [1.0, 1.0, 1.0])

% add path 
Units={'L2_left',	'L2_botm',	'L2_right',	'Icache',	'Dcache',	'Bpred',	'DTB',	'FPAdd', ...
    'FPReg',	'FPMul',	'FPMap',	'IntMap',	'IntQ',	'IntReg',	'IntExec',	'FPQ',	'LdStQ',	'ITB'};
Loc=1:18;
MTTF= mttf_em_processor(Tempr);
bar(MTTF, 'r')
axis tight;
set (gca, 'XTick', Loc);
set(gca, 'xtickLabel', Units)
set (gca, 'YLim', [0 max(MTTF)+0.1*max(MTTF)]);
xlabel('Processor Functional Units')
ylabel('MTTF in Years')

hold on 
plot(0:19, min(MTTF)*ones (1,20), 'b--', 'LineWidth', 2)

end