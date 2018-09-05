function struct2flp (FloorPlanzDb, filename, delimeter)
%given the floor plan struct or cell FloorPlanzDb, generate a flp file 
%given by filename

switch nargin
    
    case 1
        error('Incorrect inputs')
    case 2
        
        delimeter=' ';
    case 3
        
        if strcmp(delimeter, 'space')
            delm=' ';
        elseif strcmp(delimeter, '\t')
            delm=char(9);
        else
            error('Unknown delimeter')
        end
        
end

if isa(FloorPlanzDb, 'cell')
    
    label=FloorPlanzDb{:,1};
    w=FloorPlanzDb{:,2};
    h=FloorPlanzDb{:,3};
    x=FloorPlanzDb{:,4};
    y=FloorPlanzDb{:,5};
    NoOfBlocks=FloorPlanzDb{:,6};
    EmptyLabel=FloorPlanzDb{:,7};
    
    %elseif isa(FloorPlanzDb, 'struct')
    
    
else
    label=FloorPlanzDb.Labels;
    w=FloorPlanzDb.w;
    h=FloorPlanzDb.h;
    x=FloorPlanzDb.x;
    y=FloorPlanzDb.y;
    NoOfBlocks=FloorPlanzDb.NoOfBlocks;
    EmptyLabel=FloorPlanzDb.EmptyLabels;
    
    
end

%open a file with write permission
fp= fopen (filename, 'wb');
dotpoint=strfind(filename, '.');
name=filename(1:dotpoint-1);
ext=filename(dotpoint+1:end);
fid=fopen ([name, '_std.', ext], 'wb'); 

%Write the header comments
tline=['# Hotspot Floorplan file', char(10)];
fwrite(fp,tline);
fwrite(fid,tline);
tline=['# Line Format: <unit-name>\t<width>\t<height>\t<left-x>\t<bottom-y>', char(10)];
fwrite(fp,tline);
fwrite(fid,tline);
tline=['# all dimensions are in meters', char(10)];
fwrite(fp,tline);
fwrite(fid,tline);

NoOfBlocks=length(w);

for i=1:NoOfBlocks
    
    %write the data to the file
    if nargin==2
    tline=[label{i}, ' ', num2str(w(i)),' ', num2str(h(i)), ' ', ...
        num2str(x(i)),' ', num2str(y(i)), char(10)];
    else
     tline=[label{i}, delm, num2str(w(i)),delm, num2str(h(i)), delm, ...
        num2str(x(i)),delm, num2str(y(i)), char(10)];
       
    end
    format long;
    fwrite(fp,tline);
    
    fprintf(fid, '%s\t%.10f\t%.10f\t%.10f\t%.10f\n', ...
        label{i}, w(i), h(i), x(i), y(i));

    
end

fclose(fp);
fclose(fid);


end