function [FloorPlanCell FloorPlan]= flp2cell(hotspot_floor_plan)
% This function reads a standard hotspot_floor_plan and converts into a matrix
% structure for easy processing and accesing of the floor plan data


% INPUTS:
% hotspot_floor_plan = string ; file name

% OUTPUTS
% FloorPlanCell = Cell containg the fllor plan data
% FloorPlan = data stucture containg all the floor plan data

%==========================================================================
%Version    : 1.0
%Date       : July 2012
%Author     : Santanu Sarma, University of California, Irvine
%Address    : 3069, Embedded System Lab, Bren Hall, UCI, CA 92697
%Email      : santanus@uci.edu
%Web        : https://students.ics.uci.edu/~santanus/
%==========================================================================


%% Read geometry of floorplan file
% The fllor plan files can have have nay comments 

%find the name and extension of the file name
%hotspot_floor_plan='UltraSpacrT1.flp'
dotpoint=strfind(hotspot_floor_plan, '.');
name=hotspot_floor_plan(1:dotpoint-1);
ext=hotspot_floor_plan(dotpoint+1:end);


%read teh flp file
fp_src=fopen(hotspot_floor_plan, 'r')

%open a temp file without comments
newfile=[name, '.flpp']
fp= fopen (newfile, 'wb')


%get the first line
%tline = fgets(fp_src)

%check for the first char # and empty line 
while ~feof(fp_src) %ischar(tline)
     tline = fgets(fp_src);
     
%if the line is not a comment then write to the file      
    
if (tline(1)~='#')  

    %disp(tline);
    %tline = fgets(fp_src);
    
    %Write the data without comments to the new file 
    fwrite(fp, tline);

    
end

end
fclose(fp_src);
fclose(fp);

[unit_name_1, width_1, height_1, x_left_1, y_bottom_1] = ...
    textread(newfile, '%s %f %f %f %f');

NoOfBlocks=length(width_1);

% %EmptyLabels={zeros(NoOfBlocks,1)};
for i=1:NoOfBlocks
 EmptyLabels{i}='';
end

% nBlocks=length(label_mc)
% for i=1:nBlocks
% label_mc{i}='';
% end


%Compute the Width W
W=max(x_left_1 + width_1) - min(x_left_1);
%Compute the Height H
H =max(y_bottom_1 + height_1) - min(y_bottom_1);


FloorPlanCell={unit_name_1,width_1,height_1,x_left_1,y_bottom_1, ...
               NoOfBlocks,EmptyLabels, W, H};
 
%Generate a Structure of the Floorplan           
FloorPlan.BaseX=0; %coordinate of the overall core
FloorPlan.BaseY=0; %Coordinate of the overall core
           
FloorPlan.Labels= unit_name_1;
FloorPlan.w=width_1; %unit widths
FloorPlan.h=height_1; %unit heights
FloorPlan.x=x_left_1+FloorPlan.BaseX; % unit x cordinates
FloorPlan.y=y_bottom_1+FloorPlan.BaseY; %unit y cordinates
FloorPlan.NoOfBlocks=NoOfBlocks;
FloorPlan.EmptyLabels=EmptyLabels;
FloorPlan.Width=W;
FloorPlan.Height=H;

end