function solunsteady = cfdJunaidImport(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   SOLUNSTEADY = IMPORTFILE(FILENAME) Reads data from text file FILENAME
%   for the default selection.
%
%   SOLUNSTEADY = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from
%   rows STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   solunsteady = importfile('sol_unsteady.solver.stdout.dat', 3, 840);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2021/12/16 11:21:17

%% Initialize variables.
delimiter = {'   ','  ',', '};
if nargin<=2
    startRow = 3;
    endRow = inf;
end

%% Format for each line of text:
%   column2: double (%f)
%	column3: double (%f)
%   column4: double (%f)
%	column5: double (%f)
%   column6: double (%f)
%	column7: double (%f)
%   column8: double (%f)
%	column9: double (%f)
%   column10: double (%f)
%	column11: double (%f)
%   column12: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%*s%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
solunsteady = table(dataArray{1:end-1}, 'VariableNames', {'time','Rrho','Clift','Cdrag','Cdragp','Cdragv','Csidef','Cmx','Cmy','Cmz','Anglea'});

