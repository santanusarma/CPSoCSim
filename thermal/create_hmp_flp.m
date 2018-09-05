function hmp_flp = create_hmp_flp(core1_flp, core2_flp, TotalArea, NBig, NSmall)
%Not completely implemented
%Create a Heterogenous multi processor floorplan 
% core1_flp = floor plan structure of core 1 ( Big)
% core2_flp = floorplan structure of core 2 (Small)

% TotalArea= Area of the total fllor plan; Can be specified interms of the
% size of the big core

% NBig = No fo minimum big cores 
% NSmall = No of minimum small cores the 




%Test data
% single_core_flp=ev6_str_flp
% nx=2
% ny=2


unit_name_N=[];
x_left_N=zeros(nx*ny*single_core_flp.NoOfBlocks,1);
y_bottom_N=zeros(nx*ny*single_core_flp.NoOfBlocks,1);
height_N=zeros(nx*ny*single_core_flp.NoOfBlocks,1);
width_N=zeros(nx*ny*single_core_flp.NoOfBlocks,1);
l=1;
for i=1:ny %y axis
    for j=1:nx %x axis
        
        for k=1:single_core_flp.NoOfBlocks
            
            unit_name_N{l} = [char(single_core_flp.Labels(k)),...
                '_', num2str(i), '_', num2str(j)];
            EmptyLabels{l}='';
            x_left_N(l) = single_core_flp.x(k) + (j-1)*single_core_flp.Width;
            y_bottom_N(l) = single_core_flp.y(k) + (i-1)*single_core_flp.Height;
            width_N(l) = single_core_flp.w(k);
            height_N(l) = single_core_flp.h(k);
            
            l=l+1;
            
        end
        
        
    end
end

mpsoc_flp.BaseX=0;
mpsoc_flp.BaseY=0;
mpsoc_flp.Labels=unit_name_N;
mpsoc_flp.x=x_left_N;
mpsoc_flp.y=y_bottom_N;
mpsoc_flp.w=width_N;
mpsoc_flp.h=height_N;
mpsoc_flp.Width=single_core_flp.Width*nx;
mpsoc_flp.Height=single_core_flp.Height*ny;
mpsoc_flp.NoOfBlocks=single_core_flp.NoOfBlocks*nx*ny;
mpsoc_flp.EmptyLabels=EmptyLabels;


