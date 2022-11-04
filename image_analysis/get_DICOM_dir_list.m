function DICOM_folders = get_DICOM_dir_list(root_folder)
% it is understood that the file structure is root_folder/parent_folder/xxx.IMA
% this function gets a list of all DICOM folders inside of the root_folder
DICOM_folders = [];

parent_folders = dir(fullfile(root_folder));
parent_folders_bool = [parent_folders.isdir];
parent_folders = {parent_folders(parent_folders_bool).name};
parent_folders(1:2) =[]; % remove '.' and '..' from the list

for i = 1:numel(parent_folders)
    full_dir = strcat(root_folder, '/', parent_folders{i});
    DICOM_folders = [DICOM_folders, full_dir];
end
