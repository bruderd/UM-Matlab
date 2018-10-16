% snake_trysomets

t = 1 : 0.25 : 5;    % range of different t's to try
fails = [];

% load the data file to be used
[data_file,data_path] = uigetfile;
matcontents = load([data_path, data_file]); % must be a .mat file

for i = 1 : length(t)
    
    try
        [ koopman, error, data, data4sysid ] = main_snake_funoft(t(i) , matcontents);
        disp(['Worked when t = ', num2str(t(i))]);
    catch
        disp(['Failed when t = ', num2str(t(i))]);
        fails = [fails ; t(i)]; % keep track of which t values failed
    end
    
end