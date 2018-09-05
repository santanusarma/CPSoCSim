function  [scaled_flp, Scaled_flp_stuct]=scalefloorplan(floorplanfile, scalefactor)
%scale a floorplan dimensions by the scalefactor


FloorPlanzDb=flp2mat(floorplanfile);
label=FloorPlanzDb{:,1}
w=FloorPlanzDb{:,2}.*scalefactor;
h=FloorPlanzDb{:,3}.*scalefactor;
x=FloorPlanzDb{:,4}.*scalefactor;
y=FloorPlanzDb{:,5}.*scalefactor;
NoOfBlocks=FloorPlanzDb{:,6};
EmptyLabel=FloorPlanzDb{:,7};

scaled_flp= {label,w,h,x,y,NoOfBlocks,EmptyLabel};

view_floorplan(w, h, x, y, 1)

 
           
FloorPlan.Labels= label;
FloorPlan.w=w; %unit widths
FloorPlan.h=h; %unit heights
FloorPlan.x=x; % unit x cordinates
FloorPlan.y=y; %unit y cordinates
FloorPlan.NoOfBlocks=NoOfBlocks;
FloorPlan.EmptyLabels=EmptyLabel;


%Compute the Width W
W=max(x + w) - min(x);
%Compute the Height H
H =max(y + h) - min(y);


FloorPlan.Width=W;
FloorPlan.Height=H;

FloorPlan.BaseX=min(x); %coordinate of the overall core
FloorPlan.BaseY=min(y); %Coordinate of the overall core


Scaled_flp_stuct=FloorPlan;

end