function isolate_snakeInput
%isolate_snakeInput: Convert 3x1 input into 1x1 input for the snake robot
%trimmed data files
%   BE CAREFUL NOT TO USE THIS ON ANY NON-SNAKE TRIALS

% Prompt user to identify snake data file
disp(['Please select .mat file corresponding to raw data.']);
[data_file,data_path] = uigetfile;
[data_filepath,data_name,data_ext] = fileparts([data_path, data_file]);
data = load([data_path, data_file]);

% isolate input
t = data.t;
x = data.x;
u = data.u(:,2);

% save new file
fname = ['trimdataFiles', filesep, data_name, 'U', '.mat'];
save(fname, 't', 'x', 'u');

end

