function I1=view_placementmap(floorplan,temperature, Wp, Hp, hFigure, sensorloc,shape)
% floorplan = hotspot floor plan file
% temperature = steady state temperature of each block. Each column
% contains the temperature at each layers.
% Wp = No of Pixels in the Width 
% Hp = No of Pixels in the Height
% hFigure= Figure handle number to display teh thermal map
% sensorloc = a vector of sensor location indicating the block no in the
% floor plan
% shape = marker shape shpul dbe  'Plus', 'Circle', 'X-Mark'

if nargin>4
    FigNo=hFigure;
    senloc=sensorloc;
    Mark=shape
else
    FigNo=1;
    senloc=[];
    Mark='X-Mark';
end

%Read the fllorplan and its dimensions
[label, w, h, x, y] = ...
    textread(floorplan, '%s %f %f %f %f');

%% Calculate total width and height
x_left = min(x);
x_right = max(x + w);
W = x_right - x_left;
y_bottom = min(y);
y_top = max(y + h);
H = y_top - y_bottom;

NoOfBlocks= length(w); % is same as length(label);

FP=[w, h, x, y];

min_dim=min(min(FP(:,1:2)));

FPNorm=FP./min_dim;

%Step 1: define the pixle size and dimension of the image
%Wp=1024
%Hp=1024
N=Wp*Hp;  %Total No of pixel elements in the thermal map 
%Wp, Hp are no of pixles that has to map to the no 

%precision=min_dim/10;
precision=max(W,H)/max(Wp,Hp);

FPIndx=ceil(FP./precision);
wx=FPIndx(:,1);
hx=FPIndx(:,2);
xx=FPIndx(:,3);
yx=FPIndx(:,4);

%Find the center of the block to mark a sensor placement 
c_x=round(xx+wx./2)
c_y=round(yx+hx./2)

I=zeros(Wp,Hp);
for bn=1:NoOfBlocks  % for each block
    
    for i=1:wx(bn)+1  % overlap one pixel to avoid truncations issues
        for j=1: hx(bn)+1
            
        Indx_x= (xx(bn)+ i);
        Indx_y= (yx(bn)+ j);
        if (Indx_x<=Wp)&(Indx_y<=Hp)
        I(Indx_x,Indx_y)=temperature(bn,1);
        end
        end
    end
end




%Marker innsersion in the image

% markerInserter = vision.MarkerInserter('Shape','X-mark','BorderColor',...
%     'Custom','CustomBorderColor',uint8([0 0 255]));
%markerInserter = vision.MarkerInserter('Shape','X-Mark');
markerInserter = vision.MarkerInserter('Shape',Mark);

%Pts = int32([20 20; 40 40; 60 60]);
%Find the mid points of the block correpsonding to the sensorloc
Pts =[c_y,c_x];
Pts=int32(Pts);
Pts1=Pts(sensorloc,:)
%I2 = step(markerInserter, uint8(I), Pts1);
I2 = step(markerInserter, I, Pts1);
I3= imrotate(I2,90);
%find(min(I3)==0)

[B1 lims map] = real2rgb(double(I3), 'jet');
figure(FigNo)
hIm2 = imshow(B1);
set(gcf, 'Colormap', map);
set(gca, 'CLim', lims);
set(hIm2, 'CDataMapping', 'scaled');
colorbar;


%Make placement by holding the marks on the image  

I1=imrotate(I,90); %rotate the image  by 90 degree
%size(I1);
%class(I1)

[B lims map] = real2rgb(I, 'jet1');
figure(FigNo+1)
hIm = imshow(B);
set(gcf, 'Colormap', map);
set(gca, 'CLim', lims);
set(hIm, 'CDataMapping', 'scaled');
colorbar;

hold on
%plot(c_y, c_x, 'kx', 'LineWidth', 2)
plot(Pts1(:,1), Pts1(:,2), 'kx', 'LineWidth', 2)
hold off;


end