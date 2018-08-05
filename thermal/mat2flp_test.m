%test mat2flp scirpt

 clear all; clc
 
 
 % Load the EV6 fllor plan 
[unit_name_1, width_1, height_1, x_left_1, y_bottom_1] = ...
    textread('ev6_1x1.flp', '%s %f %f %f %f');
 
% Inputs:
w=width_1
h=height_1
x=x_left_1
y=y_bottom_1
hFigure=1
label=unit_name_1
z= randi(256, length(width_1),1 )

% Floor plan with out lables and colors (Pass)
view_floorplan(w, h, x, y, hFigure)



%generate a new flp file (pass)

mat2flp(w, h, x, y, label, 'ev6_1x1_new.flp')


