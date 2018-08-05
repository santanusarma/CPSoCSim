function mat2flp (w, h, x, y, label, filename)
%given the floor plan data, generate a flp file given by filename

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

%close(fp);


end