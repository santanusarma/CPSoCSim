function view_flp (floorplanfile, showlabels, figNo)
%This function visualizes teh floorplan file on figNo 
% FigNo is teh handle of teh figure to show the floorplan
% showlabels is 'showlabel'


FloorPlanzDb=flp2mat(floorplanfile);


switch (nargin)

case 1
    hFigure=1;
    label=FloorPlanzDb{:,7}; % empty lebels
   

case 2
    
    
    if strcmp(showlabels,'showlabel')
    
    label=FloorPlanzDb{:,1}

    else
        label=FloorPlanzDb{:,7}; % empty lebels
        
    end
    hFigure=1;
case 3
    hFigure=figNo;
    if strcmp(showlabels,'showlabel')
    
    label=FloorPlanzDb{:,1};

    else
        label=FloorPlanzDb{:,7}; % empty lebels
        
    end

        
    
end




w=FloorPlanzDb{:,2};
h=FloorPlanzDb{:,3};
x=FloorPlanzDb{:,4};
y=FloorPlanzDb{:,5};

view_floorplan(w, h, x, y, hFigure, label);





end