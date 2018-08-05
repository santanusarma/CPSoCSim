function cell2flp (FloorPlanzDb, filename)
%given the floor plan cell or struct FloorPlanzDb, generate a flp file 
%given by filename

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


%Write the header comments
tline=['# Hotspot Floorplan file', char(10)];
fwrite(fp,tline);
tline=['# Line Format: <unit-name>\t<width>\t<height>\t<left-x>\t<bottom-y>', char(10)];
fwrite(fp,tline);
tline=['# all dimensions are in meters', char(10)];
fwrite(fp,tline);

NoOfBlocks=length(w);

for i=1:NoOfBlocks
    
    %write the data to the file
    tline=[label{i}, ' ', num2str(w(i)),' ', num2str(h(i)), ' ', ...
        num2str(x(i)),' ', num2str(y(i)), char(10)];
    
    fwrite(fp,tline);
    
end

fclose(fp);


end