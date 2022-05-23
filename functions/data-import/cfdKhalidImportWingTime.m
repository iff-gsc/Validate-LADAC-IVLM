function [table_,time] = cfdKhalidImportWingTime(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   DISTRIBUTIONSSOL1 = IMPORTFILE(FILENAME) Reads data from text file
%   FILENAME for the default selection.
%
%   DISTRIBUTIONSSOL1 = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data
%   from rows STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   distributionssol1 = importfile('distributions_sol.surface.pval.unsteady_i=1_t=2.0000e-04.dat', 7, 1009);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2022/05/18 14:50:29


%% time
str = fileread( filename );

idx_line_end = strfind( str, newline );

str_1 = str(1:idx_line_end(1));
str_split_time = strsplit(str_1,'t=');
str_split_time_end = strsplit(str_split_time{2},'"');
time = str2double(str_split_time_end{1});

%% Initialize variables.
delimiter = {'\t',' '};
if nargin<=2
    startRow = 7;
    endRow = inf;
end

%% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%q%q%q%q%q%q%q%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
textscan(fileID, '%[^\n\r]', startRow(1)-1, 'WhiteSpace', '', 'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    textscan(fileID, '%[^\n\r]', startRow(block)-1, 'WhiteSpace', '', 'ReturnOnError', false);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,2,3,4,5,6,7]
    % Converts text in the input cell array to numbers. Replaced non-numeric
    % text with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^[-/+]*\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end


%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
not_nan_idx = find(~isnan(cell2mat(raw(:,1))));
table_ = table;
table_.Y = cell2mat(raw(not_nan_idx, 1));
table_.fx = cell2mat(raw(not_nan_idx, 2));
table_.fy = cell2mat(raw(not_nan_idx, 3));
table_.fz = cell2mat(raw(not_nan_idx, 4));
table_.mx = cell2mat(raw(not_nan_idx, 5));
table_.my = cell2mat(raw(not_nan_idx, 6));
table_.mz = cell2mat(raw(not_nan_idx, 7));
table_.eta = table_.Y/max(table_.Y);


end
