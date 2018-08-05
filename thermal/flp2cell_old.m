function FloorPlanCell= flp2cell_old(hotspot_floor_plan)
% This function reads a hotspot_floor_plan and converts into a cell
% structure for easy processing and accesing of the floor plan data


% INPUTS:
% hotspot_floor_plan = string ; file name

% OUTPUTS
% FloorPlanCell = Cell containg the fllor plan data

%==========================================================================
%Version    : 1.0
%Date       : July 2012
%Author     : Santanu Sarma, University of California, Irvine
%Address    : 3069, Embedded System Lab, Bren Hall, UCI, CA 92697
%Email      : santanus@uci.edu
%Web        : https://students.ics.uci.edu/~santanus/
%==========================================================================


%% Read geometry of floorplan file
% The fllor plan files should not have nay comments 

[unit_name_1, width_1, height_1, x_left_1, y_bottom_1] = ...
    textread(hotspot_floor_plan, '%s %f %f %f %f');

FloorPlanCell={unit_name_1,width_1,height_1,x_left_1,y_bottom_1};
end