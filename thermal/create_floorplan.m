% CREATE_FLOORPLAN creates a new multi-core floorplan by tiling copies of 
% a single core floorplan.
%
%   CREATE_FLOORPLAN(single_core_FP_file, N_x, N_y, multi_core_FP_file) 
%   creates a multi-core floorplan with N_x copies of the the single core 
%   floorplan (specified in the file single_core_FP_file) horizontally and
%   N_y copies vertically. This floorplan is then scaled to fit in the same
%   dimensions as the original single core floorplan. (The file format for 
%   the floorplan file is explained in the comments at the beginning of the
%   example floorplan file ev6.flp.bak.) The resulting multi-core floorplan 
%   is stored in the file multi_core_FP_file.
%
% Ravishankar Rao, Arizona State University
% Created, Dec 6, 2006

% INPUTS:
% single_core_FP_file = string, containing hotspot floor plan file with *.flp 
% N_x = No of cores in x direction
% N_y = No of cores in y direction
% multi_core_FP_file = String; Name of the multicore hotspot fllr plan with *.flp
% extesnion 

% OUTPUTS:
% generate the file multi_core_FP_file with the floor plan values. Fo rthe
% format check hotspot floor plan file 

% EXAMPLE:
% create_floorplan ('ev6_1x1.flp', 2,2, 'mc_ev6_2x2.flp')
%
% Will cerate a 4 core (2x2) floor plan in the file 'mc_ev6_2x2.flp' form
% the singel core fllor plan 'ev6_1x1.flp'.
 

function create_floorplan(single_core_FP_file, N_x, N_y, multi_core_FP_file)
    
%% Read geometry of single core floorplan
[unit_name_1, width_1, height_1, x_left_1, y_bottom_1] = ...
    textread(single_core_FP_file, '%s %f %f %f %f');

%% Compute chip width and height
core_width_1  = max(x_left_1 + width_1) - min(x_left_1);
core_height_1 = max(y_bottom_1 + height_1) - min(y_bottom_1);
N_units_per_core = length(x_left_1);
chip_width = core_width_1;
chip_height = core_height_1;

%% Compute dimensions of each core of multicore floorplan
core_width = chip_width/N_x;
core_height = chip_height/N_y;
scale_factor_x = 1/N_x;
scale_factor_y = 1/N_y;
width = width_1*scale_factor_x;
height = height_1*scale_factor_y; 

%% Compute locations of functional units for bottom left core
x_left = x_left_1*scale_factor_x;
y_bottom = y_bottom_1*scale_factor_y + (N_y-1)*core_height;  

%% Create multicore floorplan file
fid = fopen(multi_core_FP_file, 'w');  

%% Compute locations of each core as shifted versions of the bottom left
% core and write this information to file in .flp file format.
for i = 1:N_y
    for j = 1:N_x
        % For each functional unit k in core (i, j)
        for k = 1:N_units_per_core
            unit_name_N = char(unit_name_1(k));
            x_left_N = x_left(k) + (j-1)*core_width;
            y_bottom_N = y_bottom(k) - (i-1)*core_height;
            width_N = width(k);
            height_N = height(k);
            fprintf(fid, '%s_%d_%d\t%.10f\t%.10f\t%.10f\t%.10f\n', ...
                unit_name_N, i, j, width_N, height_N, x_left_N, ...
                y_bottom_N);
        end
    end
end  

%% Close file
fclose(fid);
