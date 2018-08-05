function [joined_ptrace joined_flp]=join_ptrace(ptrace1, flp1, ptrace2, flp2, joined_ptrace_file, position)
%This function joins two ptrace files or structures for  the given
%position based on the floorplan and returns the joined ptrace file and
%joined flp file


% Eg:
% flp1= 'mpsoc_ev4_3x3.flp' 
% flp2= 'ev6_hetero.flp'
% ptrace1 = './power_data/ptrace/mpsoc_ev4_3x3.ptrace'
% ptrace2 = './power_data/ptrace/gcc.ptrace'
% joined_ptrace_file='joined.ptarce'
% position='left'

%ptrace1
%ptrace2

switch (nargin)
    case 1
        error('invalid input arguments');
        
    case 2
        error('invalid input arguments');
        
    case 3 % ptrace1, ptrace2, joined_ptrace_file
        
        %Get the labels from the ptrace file
        if isa(ptrace1, 'char')
            power_tr1=ptrace2pstruct(ptrace1)
            labels1=power_tr1.textdata   %labels
        else
             error('invalid ptrace1 format');
        end
        
        if isa(ptrace2, 'char')
            power_tr2=ptrace2pstruct(ptrace2);
            labels2=power_tr2.textdata  ; %labels
        else
            error('invalid ptrace2 format');
        end

        position='left'  
        NewLabels ={labels1{:}, labels2{:}}
        
    case 4
        error('invalid input arguments');
    case 5
        position='left'
    case 6
        % Read flp 1 and flp2 to get the labels
        if isa ( flp1, 'char')
            [~ , flp1_str] =flp2struct(flp1);
        else
            error('invalid flp1 format');
        end
        
        if isa ( flp2, 'char')
            [~ , flp2_str] =flp2struct(flp2);
        else
            error('invalid flp2 format');
        end
        
        
        %Get the labels from the ptrace file
        if isa(ptrace1, 'char')
            power_tr1=ptrace2pstruct(ptrace1)
            labels1=power_tr1.textdata   %labels
        else
            error('invalid ptrace1 format');
        end
        
        if isa(ptrace2, 'char')
            power_tr2=ptrace2pstruct(ptrace2)
            labels2=power_tr2.textdata   %labels
        else
            error('invalid ptrace2 format');
        end
        
        %NewLabels ={labels1{:}, labels2{:}}
        %NewLabels= {flp1_str.Labels{:}, flp2_str.Labels{:}};
        
        %joined_flp = join_flp(flp1_str, flp2_str, position);
        
        joined_flp = join_flps(flp1_str, flp2_str, position);
        
        %Adjust the base
        
        
        NewLabels= joined_flp.Labels;
        
    otherwise 
        error('invalid input arguments');
    
end
    






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

%Check the dimension of the power traces before combining

power= [power_tr1.data(1:NoOfSamplingPoints,:),power_tr2.data(1:NoOfSamplingPoints,:)];

joined_ptrace.data=power;
for i=1:NoOfSamplingPoints
    
    tline=[num2str(power(i,:)), char(10)];
    fwrite(fp1,tline);
end

fclose(fp1);

end