function vect = mat2vect (A, dim_rows, dim_cols)
%This function convects a (2D) matrix A into a vector (1D) vect


if nargin==1
    [rows, cols]= size(A);
    
    siz=rows*cols;
    vect=reshape(A, siz, 1);
    
elseif nargin==2 %check this case 
    rows=dim_rows;
    cols=dim_rows;
    siz=rows*cols;
    vect=reshape(A, siz, 1);
        
elseif nargin==3
    
    rows=dim_rows;
    cols=dim_cols;
    [sx, sy]=size(A);
    siz=rows*cols;
        
        if(siz==sx*sy)
        
        vect=reshape(A, siz, 1);
        else
            error('Incorrect dimension');
        end
       
else
    error('Unknown format');
end


end