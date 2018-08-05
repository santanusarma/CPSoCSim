function view_grid_ttrace(grid_ttrace_file, dim_rows, dim_cols, plottype, units, varargin)
%This function views the thermal trace of the grid model output of Hotspot
%check the x, y directions ( node maddping in 2D) and map the floorplan

%The default trace values are in ^oC

if nargin==1
    x=vect2mat(grid_ttrace_file);
    surf(x);
    temp_unit='C';
    plottype='surf';
elseif nargin==3
    temp_unit='C';
    plottype='surf';
    x=vect2mat(grid_ttrace_file, dim_rows, dim_cols);
    surf(x);
elseif nargin==4
    temp_unit='C';
    x=vect2mat(grid_ttrace_file, dim_rows, dim_cols);
    switch plottype
        case 'surf'
            surf(x);
            
        case 'mesh'
            mesh(x);
            
        case 'waterfall'
            waterfall(x);
            
        otherwise
            surf(x);
            
            
    end
    xlabel('x index')
    ylabel('y index')
    
elseif nargin==5
    
    switch units
        
        case 'Centigrade'
            temp_unit='C';
            x=vect2mat(grid_ttrace_file, dim_rows, dim_cols)...
                -273.15*ones(dim_rows,dim_cols);
        case 'Kelvin'
            temp_unit='K';
            x=vect2mat(grid_ttrace_file);
        otherwise
            temp_unit='C';
            x=vect2mat(grid_ttrace_file, dim_rows, dim_cols);
            
    end
    switch plottype
        case 'surf'
            surf(x);
            
        case 'mesh'
            mesh(x);
            
        case 'waterfall'
            waterfall(x);
            
        otherwise
            surf(x);
            
            
    end
elseif nargin==6
    
    figno=varargin{1}
    figure(figno)
    switch units
        
        case 'Centigrade'
            temp_unit='C';
            x=vect2mat(grid_ttrace_file, dim_rows, dim_cols)...
                -273.15*ones(dim_rows,dim_cols);
        case 'Kelvin'
            temp_unit='K';
            x=vect2mat(grid_ttrace_file);
        otherwise
            temp_unit='C';
            x=vect2mat(grid_ttrace_file, dim_rows, dim_cols);
            
    end
    switch plottype
        case 'surf'
            surf(x);
            
        case 'mesh'
            mesh(x);
            
        case 'waterfall'
            waterfall(x);
            
        otherwise
            surf(x);
            
            
    end
    
    
else
    error('Incorrect Format')
    
end

FontSize=16;
zlabel(['Temperature, ^o',  temp_unit], 'FontSize', FontSize, 'FontWeight', 'bold','Color', 'k')
xlabel('x-index', 'FontSize', FontSize, 'FontWeight', 'bold','Color', 'k')
ylabel('y-index', 'FontSize', FontSize, 'FontWeight', 'bold','Color', 'k')
%axis([1 64 1 64 70 110])
set (gcf, 'Color', [1.0, 1.0, 1.0])
set(gca, 'LineWidth', 2)
set(gca, 'FontSize', FontSize-2)
set(gca, 'FontWeight', 'bold')



end