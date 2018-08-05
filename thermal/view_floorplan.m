% VIEW_FLOORPLAN provides a visual representation of a floorplan with
% optional labels for each block and a 'z-value' for each block that is
% displayed as varying levels of grayscale color. 
%
%   VIEW_FLOORPLAN(w, h, x, y, hFigure) plots a floorplan with N 
%   blocks, withthe i^th block having width w_i, height h_i, having its 
%   bottom left corner at (x_i, y_i). w is a column vector consisting of  
%   all block widths i = 1, ..., N. Similarly, for vectors x, y, and h. 
%   hFigure is the handle to the figure object used for the plot.
%
%   VIEW_FLOORPLAN(w, h, x, y, hFigure, label) adds a text label for
%   each block, where the labels are specified as a cell array of strings.
%
%   VIEW_FLOORPLAN(w, h, x, y, hFigure, label, z) shows the z-value
%   for each block as a grayscale level, with black representing the least
%   value and white the highest.
%
%   See test_view_floorplan.m for an example on using this function.
%
% Ravishankar Rao, Arizona State University
% Created, Oct 12 2007

%==========================================================================
% ADDITIONAL HELP & TEST CASES
%==========================================================================
% Version    : 1.0
% Date       : July 2012
% Author     : Santanu Sarma, University of California, Irvine
% Address    : 3069, Embedded System Lab, Bren Hall, UCI, CA 92697
% Email      : santanus@uci.edu
% Web        : https://students.ics.uci.edu/~santanus/
%==========================================================================
% INPUTS :
% w = width of each block, an 1d array 
% h = hight of each block, an 1d array 
% x = bottom left corner, x cordinate, an 1d array
% y = bottom left corner, y cordinate, an 1d array
% z = a 1d array of bw color (in 0 to 256 range)  

% OUTPUTS :
% A figure showing the floor plan given by the figure handle hFigure

% EXAMPLE:
% See view_floorplan_test.m
%==========================================================================

function view_floorplan(w, h, x, y, hFigure, label, z)

%% Show z value through color
if nargin > 6
    showZ = 1;
    z_min = min(z);
    z_max = max(z);
    z_norm = (z - z_min)/(z_max - z_min);
    % Grayscale: z_min = black, z_max = white
    z_color = repmat(z_norm, 1, 3); 
else
    showZ = 0;
end

if (nargin > 5) && ~isempty(label) 
    showLabel = 1; 
else
    showLabel = 0; 
end

%% Setup figure
figure(hFigure);
set (hFigure, 'Color', [1.0, 1.0, 1.0])
clf;
axis equal;
set(gca, 'FontSize', 16);

%% Calculate total width and height
x_left = min(x);
x_right = max(x + w);
W = x_right - x_left;
y_bottom = min(y);
y_top = max(y + h);
H = y_top - y_bottom;

%% Draw bounding rectangle
rectangle('Position', [x_left y_bottom W H], 'EdgeColor', 'blue', 'LineWidth', 2);

%% Draw rectangles for individual blocks
N_blocks = length(x);
for i = 1:N_blocks
    if showZ == 1
        rectangle('Position', [x(i) y(i) w(i) h(i)], 'FaceColor', ...
            z_color(i,:), 'EdgeColor', 'blue','LineWidth', 2);
    else
        rectangle('Position', [x(i) y(i) w(i) h(i)], 'EdgeColor', 'blue','LineWidth', 2);
    end
    if showLabel== 1
        text(x(i)+w(i)/3, y(i)+h(i)/3, label{i}, 'BackgroundColor',...
            [.7 .9 .7], 'FontSize', 16,'LineWidth', 2);
    end
    
    %pause(2)
end
axis tight;

%% Export plot to PDF
print -dpdf floor_plan_figure;

