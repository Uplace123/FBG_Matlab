%% Note to self:
% "20221028": PSM, l = 20, 40, 60

%%
addpath(genpath('helper_funcs'));

DICOM_root_folder = "~/Desktop/DICOM/20221104"; % Where the root folder is located
data_specs_folder = "~/Documents/FBG_Matlab/image_analysis";
data_specs_filename = 'data_specs.csv'; % output file
threshold = 15000; % intensity-based threshold for segmentation


data_specs_filepath = fullfile(data_specs_folder, data_specs_filename);
if exist(data_specs_filepath, 'file') == 2
    warning('Previous data spec sheet detected.')
%     delete(data_specs_filepath);
return;
end


DICOM_folders = get_DICOM_dir_list(DICOM_root_folder);
DICOM_folders = string(DICOM_folders)'; 

skip_datasets = []; % problematic datasets to be skipped
if ~isempty(skip_datasets)
    skip_folders = strcat(DICOM_root_folder, '/', skip_datasets);
    for i = 1:numel(skip_folders)
        warning(sprintf("Skipping dataset %s\n", skip_folders(i)));
    end
    read_folders = setdiff(DICOM_folders, skip_folders, 'stable');
else
    skip_folders = [];
    read_folders = DICOM_folders;
end

sbs = zeros(size(read_folders));
ls = zeros(size(read_folders));
dbs = zeros(size(read_folders));

tot_folders = numel(read_folders);
for i = 1:tot_folders
    fprintf('(%d/%d) Reading %s: ', i, tot_folders, read_folders(i));
    l = group_folders_func(read_folders(i), threshold);
    ls(i) = round(l/5)*5; % round to nearest 5
    fprintf('l = %d\n', ls(i));
end

output_table = table(read_folders, ls);
writetable(output_table, data_specs_filepath)
disp("FINISHED")



%% group_folders_func
function l = group_folders_func(dataset, threshold)
% Read volume
small_vol_threshold = 40;
[newVol, delMM] = readCT(dataset, threshold, small_vol_threshold);

% Read frames
% Gel Frame
gel_marker_design_coord = readMarkerDesignCoord("gel_marker_coord.txt");
n_gel_markers = size(gel_marker_design_coord, 1);

% Assign ID
CC = bwconncomp(newVol); % find connected regions
numPixels = cellfun(@numel,CC.PixelIdxList);
[~, needle_id] = maxk(numPixels, 1); % needle will correspond to the largest group of connected voxels
[~, gel_markers_id] = maxk(numPixels, 1 + n_gel_markers);
gel_markers_id = setdiff(gel_markers_id, needle_id); % gel markers will correspond to next n_gel_markers largetst groups of voxels
totList = linspace(1, CC.NumObjects, CC.NumObjects);
stage_markers_id = setdiff(totList, [needle_id, gel_markers_id]); % stage markers will be the last few groups; it is set last because some markers might be invisible

% Segmentation
% Voxel to MM scale
V2MM_scale = delMM*eye(3);

% Gel markers
gelMarkerVol = segmentVolume(size(newVol), CC, gel_markers_id);
gelMarkerCentroidsTab = regionprops3(gelMarkerVol, 'Centroid');
gel_marker_CT_coord = table2array(gelMarkerCentroidsTab)*V2MM_scale;
% volumeViewer(gelMarkerVol); % Visualize the volume

% Needle
needleVol = segmentVolume(size(newVol), CC, needle_id);
needleVolSkel = bwskel(needleVol); % use the skeleton approach here
[needleVoxI, needleVoxJ, needleVoxK] = findND(needleVolSkel);
needleMM = [delMM*needleVoxJ, delMM*needleVoxI, delMM*needleVoxK];% For some reason the I and J are flipped


% Registration
% Gel -> Find (R_gel, t_gel) to bring volume to design
[R_gel, t_gel, err_gel] = icp_wrapper(gel_marker_design_coord, gel_marker_CT_coord, 1);
gel_marker_CT_coord = apply_transform(gel_marker_CT_coord, R_gel, t_gel);

% Needle
needleMM = apply_transform(needleMM, R_gel, t_gel);

% Output
l = max(needleMM(:, 1));
end
