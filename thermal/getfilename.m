function [filename, ext] =getfilename(file_path)
% Get the file name and the extension hen the input string contain [path
% and the file name

% eg flie_path= '../../test.txt'
% filename= test
% ext=.txt

pathstr=strfind(file_path, '/');
dotpoint=strfind(file_path, '.');
if isempty(pathstr) && length(dotpoint)==1
    filename=file_path(1:dotpoint-1);
    ext=file_path(dotpoint+1:end);
else
    filename=file_path(pathstr(end)+1:dotpoint(end)-1);
    ext=file_path(dotpoint(end)+1:end);
    
end