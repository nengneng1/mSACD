% This function calculates the list of folder trees in an address 
function folderAddresses = findEndSubfolders(parentFolder)
    % Initialize an empty cell array to store folder addresses
    folderAddresses = {};

    % Get a list of all items within the parent folder
    items = dir(parentFolder);

    % Iterate through each item
    for i = 1:numel(items)
        % Ignore current and parent directories
        if strcmp(items(i).name, '.') || strcmp(items(i).name, '..')
            continue;
        end

        % Check if the item is a folder
        if items(i).isdir
            % Get the full path of the subfolder
            subfolderPath = (fullfile(parentFolder, items(i).name));

            % Recursively call the function for subfolders
            subfolderAddresses = findEndSubfolders(subfolderPath);

            % If there are no further subfolders, add the current subfolder to the addresses
            if isempty(subfolderAddresses)
                folderAddresses{end+1} = (subfolderPath);
            else
                % Append the subfolder addresses to the existing addresses
                folderAddresses = [folderAddresses (subfolderAddresses)];
            end
        end
    end
end