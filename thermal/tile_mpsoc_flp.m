function mpsoc_flp  = tile_mpsoc_flp(single_core_flp, nx, ny)
%Create a homogeneous MPSOC floorplan 
%single_corre_flp = structure based floorplan of the single core floorplan
%that follows the same output as flp2cell
%nx = no of cores in the x axis
%ny= no of cores in the uy axis

% mpsoc_flp = structuree containing the floorplan data
% mpsoc_flp_file = file name of the mpsoc floorplan

%The output floorppan is is also in the structure format

%Should also update for file based fllorplan 



%Test data
% single_core_flp=ev6_str_flp
% nx=2
% ny=2

if isa(single_core_flp, 'char')
    [~, single_core_flp]=flp2struct(single_core_flp);
end

unit_name_N=[];
x_left_N=zeros(nx*ny*single_core_flp.NoOfBlocks,1);
y_bottom_N=zeros(nx*ny*single_core_flp.NoOfBlocks,1);
height_N=zeros(nx*ny*single_core_flp.NoOfBlocks,1);
width_N=zeros(nx*ny*single_core_flp.NoOfBlocks,1);
l=1;
for i=1:nx %y axis
    for j=1:ny %x axis
        
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
mpsoc_flp.x=x_left_N+mpsoc_flp.BaseX;
mpsoc_flp.y=y_bottom_N+mpsoc_flp.BaseY;
mpsoc_flp.w=width_N;
mpsoc_flp.h=height_N;
mpsoc_flp.Width=single_core_flp.Width*nx;
mpsoc_flp.Height=single_core_flp.Height*ny;
mpsoc_flp.NoOfBlocks=single_core_flp.NoOfBlocks*nx*ny;
mpsoc_flp.EmptyLabels=EmptyLabels;

% mpsoc_flp.nx=nx;
% mpsoc_flp.ny=ny;
% mpsoc_flp.TotalCores=nx*ny;

