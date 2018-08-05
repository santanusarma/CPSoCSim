function [FloorPlanCell FloorPlan]= flp2struct(hotspot_floor_plan,shift,Cord, scale, sf)
% This function reads a standard hotspot_floor_plan and converts into a flp
% structure for easy processing and accesing of the floor plan data


% INPUTS:
% hotspot_floor_plan = string ; file name with extension; 
% shift is string 'shift' that translates the origin by Cord [BaseX BaseY]
% Cord = a vector [ BaseX BaseY] that shift the origin
% scale is a string 'scale' which scales the flp by the scalefactor sf wrt
% to centre 
% sf = scalefactor, +ve number or fraction used for scaling the flp

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
switch nargin 
    case 1
        BaseX=0;
        BaseY=0;
        sf=1.0;
    case 2
        BaseX=0;
        BaseY=0;
        sf=1.0;
      
    case 3
        if strcmp(shift,'shift')
        BaseX=Cord(1);
        BaseY=Cord(2);
        sf=1.0;
        else
            error('error in arguments')
            
        end

    case 4 % check 
        if strcmp(shift,'shift')
        BaseX=Cord(1);
        BaseY=Cord(2);
        sf=1.0;
        elseif strcmp(shift,'scale')
        BaseX=0;
        BaseY=0;
        sf=sf;
            
        else
            error('error in arguments')
        end
        
    case 5
        if strcmp(shift,'shift') && strcmp(scale,'scale')
        BaseX=Cord(1)
        BaseY=Cord(2)
        sf=sf;
            
        else
            error('unknown arguments')
        end

end
    
dotpoint=strfind(hotspot_floor_plan, '.');
name=hotspot_floor_plan(1:dotpoint-1);
ext=hotspot_floor_plan(dotpoint+1:end);


%read teh flp file
fp_src=fopen(hotspot_floor_plan, 'r');

%open a temp file without comments
newfile=[name, '.flpp']
fp= fopen (newfile, 'wb');


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

EmptyLabels={zeros(NoOfBlocks,1)};
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


% FloorPlanCell={unit_name_1,width_1,height_1,x_left_1,y_bottom_1, ...
%                NoOfBlocks,EmptyLabels, W, H};

%Generate a Structure of the Floorplan  

% FloorPlan.BaseX=BaseX; %coordinate of the overall core
% FloorPlan.BaseY=BaseY; %Coordinate of the overall core
%Modified on July 12 to correct bug in cpsoc floorplan
FloorPlan.BaseX=BaseX+min(x_left_1); %coordinate of the overall core
FloorPlan.BaseY=BaseY+min(y_bottom_1); %Coordinate of the overall core
 
           
FloorPlan.Labels= unit_name_1;
FloorPlan.w=width_1.*sf; %unit widths
FloorPlan.h=height_1.*sf; %unit heights
FloorPlan.x=x_left_1.*sf+FloorPlan.BaseX; % unit x cordinates
FloorPlan.y=y_bottom_1.*sf+FloorPlan.BaseY; %unit y cordinates
FloorPlan.NoOfBlocks=NoOfBlocks;
FloorPlan.EmptyLabels=EmptyLabels;
FloorPlan.Width=W.*sf;
FloorPlan.Height=H.*sf;



FloorPlanCell={unit_name_1,width_1.*sf,height_1.*sf,...
                x_left_1.*sf+FloorPlan.BaseX,...
                y_bottom_1.*sf+FloorPlan.BaseY, ...
               NoOfBlocks,EmptyLabels, W.*sf, H.*sf,BaseX,BaseY};


end