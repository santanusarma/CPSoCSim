function I1=view_thermalmap(floorplan,temperature, Wp, Hp, hFigure)
% floorplan = hotspot floor plan file
% temperature = steady state temperature of each block. Each column
% contains the temperature at each layers.
% Wp = No of Pixels in the Width 
% Hp = No of Pixels in the Height
% hFigure= Figure handle number to display teh thermal map

if nargin>4
    FigNo=hFigure;
else
    FigNo=1;
end
[label, w, h, x, y] = ...
    textread(floorplan, '%s %f %f %f %f');

%% Calculate total width and height
x_left = min(x);
x_right = max(x + w);
W = x_right - x_left;
y_bottom = min(y);
y_top = max(y + h);
H = y_top - y_bottom;


FP=[w, h, x, y];

min_dim=min(min(FP(:,1:2)));

FPNorm=FP./min_dim;

%Step 1: define the pixle size and dimension of the image
%Wp=1024
%Hp=1024
N=Wp*Hp;  %No of elements in the thermal map 
%Wp, Hp are no of pixles that has to map to the no 

%precision=min_dim/10;
precision=max(W,H)/max(Wp,Hp);

FPIndx=ceil(FP./precision);
wx=FPIndx(:,1);
hx=FPIndx(:,2);
xx=FPIndx(:,3);
yx=FPIndx(:,4);

I=zeros(Wp,Hp);
for bn=1:30  % for each block
    
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

I1=imrotate(I,90); %rotate the image  by 90 degree
size(I1)
[B lims map] = real2rgb(I1, 'jet1');
figure(FigNo)
hIm = imshow(B);
set(gcf, 'Colormap', map);
set(gca, 'CLim', lims);
set(hIm, 'CDataMapping', 'scaled');
colorbar;


end