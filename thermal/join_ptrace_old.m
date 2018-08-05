function joined_ptrace =join_ptrace(ptrace1, flp1, ptrace2, flp2, joined_ptrace_file)


% Eg:
% flp1= 'ev6_hetero.flp'
% flp2= 'ev6_hetero.flp'
% ptrace1 = './power_data/ptrace/gcc.ptrace'
% ptrace2 = './power_data/ptrace/gcc.ptrace'
% joined_ptrace_file='joined.ptarce'

% Read flp 1  


if isa ( flp1, 'char')
    [~ , flp1_str] =flp2struct(flp1);
end

if isa ( flp2, 'char')
    [~ , flp2_str] =flp2struct(flp2);
end


NewLabels= {flp1_str.Labels{:}, flp2_str.Labels{:}};
joined_ptrace.Labels=NewLabels;

%Read the power traces
power_tr1=ptrace2pstruct(ptrace1)
power_tr2=ptrace2pstruct(ptrace2)

[px1, py1]=size(power_tr1.data)
[px2, py2]=size(power_tr2.data)
NoOfBlocks = py1 + py2;

%Generate trace file for the dynamic power for specified sampling points
filename=joined_ptrace_file;
fp1= fopen (filename, 'w')

%write the unit names in the first line

for i=1:NoOfBlocks
    tline=[NewLabels{i}, char(9)];
    fwrite(fp1,tline);
end
fwrite(fp1, char(10));

%compute the no of sampling point

NoOfSamplingPoints= min(px1,px2)
power= [power_tr1.data,power_tr2.data];
joined_ptrace.data=power;
for i=1:NoOfSamplingPoints
    
    tline=[num2str(power(i,:)), char(10)];
    fwrite(fp1,tline);
end

fclose(fp1);

end