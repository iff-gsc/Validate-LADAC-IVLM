function [table_3d,cell_of_tables] = cfdJunaidImportWingAll(folder_name)

disp('RANS import started.')
tic

% Find folder name in path
matlab_path = path;
len_name = length(folder_name);
idx_on_path = strfind( matlab_path, folder_name );
dir_len = 1;
for i = 1:length(idx_on_path)
    current_str = matlab_path(idx_on_path(i) + len_name);
    if ~( strcmp(current_str,';') || strcmp(current_str,':') )
        if i == length(idx_on_path)
            error('folder not detected')
        end
        continue;
    end
    while true
        current_str = matlab_path(idx_on_path(i) - dir_len);
        if ( strcmp(current_str,';') || strcmp(current_str,':') ) || idx_on_path(i) - dir_len <= 1
            if idx_on_path(i) - dir_len > 1
                dir_len = dir_len - 1;
            end
            break;
        else
            dir_len = dir_len + 1;
        end
    end
    break;
end
idx_begin = idx_on_path(i) - dir_len;
idx_end = idx_on_path(i) + len_name - 1;
% This is the path of the specified folder name
folder_path = matlab_path( idx_begin:idx_end );

% Find all files inside folder
dat_files = dir(folder_path);

num_files = length(dat_files);

cell_of_tables = {};
table_3d_T = table;
time_vec = [];

for i = 1:num_files
    if contains(dat_files(i).name,'.dat') && contains(dat_files(i).name,'i=')
        [table_,time] = cfdJunaidImportWingTime(dat_files(i).name);
        var_names = table_.Properties.VariableNames;
        if ~isempty(table_3d_T)
            depth = size(table_3d_T.(var_names{1}),2);
            for j = 1:length(var_names)
                table_3d_T.(var_names{j})(:,depth+1) = table_.(var_names{j});
            end
        else
            table_3d_T = table_;
        end
        time_vec(end+1) = time;
        
        cell_of_tables{i,1} = table_;
        cell_of_tables{i,2} = time;

        % disp(['t=',num2str(time),'s']);
        
    end
end

[time_vec_sort,idx_sort] = sort(time_vec);

% sort and transpose
table_3d = table;
for i = 1:length(var_names)
	table_3d.(var_names{i}) = table_3d_T.(var_names{i})(:,idx_sort)';
end
table_3d.time = time_vec_sort';


elapsed_time = toc;
disp(['RANS import finished after ',num2str(elapsed_time),' seconds.'])


end