function A = vect2mat(in_vec_file, grid_rows, grid_cols)
%Given a input vector or file of the grid model, convert to equivalent 
%2D matrix 

%The fault processing is for the the Hotspot v4.0 

if nargin==1
    if isa(in_vec_file, 'numeric')
        grid_data=in_vec_file;
         vector=grid_data(:,1);
    elseif isa(in_vec_file, 'char')
        grid_data=load(in_vec_file);
        
        vector=grid_data(:,2);
    else
      error('Incorrect Format')  
        
    end
      
    rows=50;
    cols=50;

elseif nargin==2
    error('Incorrect Format')

elseif nargin==3
    
    if isa(in_vec_file, 'numeric')
        grid_data=in_vec_file;
        vector=grid_data(:,1);
    elseif isa(in_vec_file, 'char')
        grid_data=load(in_vec_file);
        
        vector=grid_data(:,2);
       
    else
      error('Incorrect Format')  
        
    end
    rows=grid_rows;
    cols=grid_cols;
    
else
    error('Unknown Fromat')
end



     A=reshape(vector,rows,cols);

end