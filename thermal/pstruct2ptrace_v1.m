function pstruct2ptrace_v1 (mpsoc_pstruct, filename)
%genrate a ptrace file from pstruct power data structure 
%mpsoc_pstruct = power trace structure format
%filename = filename to write the power trace

% dotpoint=strfind(filename, '.');
% name=filename(1:dotpoint-1);
% ext=filename(dotpoint+1:end);
% newfilename=[name, '.ptrace']


%Generate trace file for the dynamic power for specified sampling points
fp1= fopen (filename, 'w')
%write the unit names in the first line

[px py]= size(mpsoc_pstruct.data)
NoOfBlocks = py;
NoOfSamplingPoints=px;

for i=1:mpsoc_pstruct.NoOfBlocks
    tline=[mpsoc_pstruct.Labels{i}, char(9)];
    fwrite(fp1,tline);
end
fwrite(fp1, char(10));

%NoOfSamplingPoints=mpsoc_pstruct.NoOfSamplingPoints;

for i=1:NoOfSamplingPoints
    
    tline=[num2str(mpsoc_pstruct.power(i,:)), char(10)];
    fwrite(fp1,tline);
end

fclose(fp1);
end