function saveFigureAsPNG(resolution, filename)
    % Check if the 'Figures' folder exists, create it if not
    workingDirectory = pwd;  % Get the current working directory
    folderName = 'Figures';
    directory = fullfile(workingDirectory, folderName);  % Create the full path to the 'Figures' folder
    if ~exist(directory, 'dir')
        mkdir(directory);
    end
    % Save the current figure as a PNG
    fullFileName = fullfile(directory, filename);
    exportgraphics(gcf, fullFileName, 'Resolution', resolution);
end