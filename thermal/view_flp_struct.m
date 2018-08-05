function view_flp_struct(floorplan_struct, showlabels, figNo)
%This function visualizes teh floorplan file on figNo 
% FigNo is the handle of the figure to show the floorplan
% showlabels is 'showlabel' to show the labels in te figure


%FloorPlanzDb=flp2mat(floorplanfile);


switch (nargin)

case 1
    hFigure=1;
    label=floorplan_struct.EmptyLabels; % empty lebels
   

case 2
    
    
    if strcmp(showlabels,'showlabel')
    
    label=floorplan_struct.Labels;

    else
        label=floorplan_struct.EmptyLabels; % empty lebels
        
    end
    hFigure=1;
case 3
    hFigure=figNo;
    if strcmp(showlabels,'showlabel')
    
    label=floorplan_struct.Labels;

    else
        label=floorplan_struct.EmptyLabels; % empty lebels
        
    end

        
    
end




w=floorplan_struct.w;
h=floorplan_struct.h;
x=floorplan_struct.x;
y=floorplan_struct.y;

view_floorplan(w, h, x, y, hFigure, label);





end