function [Nel, t_stage, y_fit_obj] = analyzeBaselineImage(dataset, threshold)
% Yanzhou Wang
% Jan 2022

%% Add directory
path = fullfile(matlabdrive, 'needle_tissue_FEM_model/image_analysis');
addpath(genpath(path));

%% Read CT
small_vol_threshold = 40;
[newVol, delMM] = readCT(dataset, threshold, small_vol_threshold);
% volumeViewer(newVol);

%% Read Body Frames
% Stage Frame
stage_marker_design_coord = readMarkerDesignCoord("stage_marker_coord.txt");
n_stage_markers = size(stage_marker_design_coord, 1);

% Gel Frame
gel_marker_design_coord = readMarkerDesignCoord("gel_marker_coord.txt");
n_gel_markers = size(gel_marker_design_coord, 1);

%% Assign ID based on connected voxels
CC = bwconncomp(newVol); % find connected regions
numPixels = cellfun(@numel,CC.PixelIdxList);
[~, needle_id] = maxk(numPixels, 1); % needle will correspond to the largest group of connected voxels
[~, gel_markers_id] = maxk(numPixels, 1 + n_gel_markers);
gel_markers_id = setdiff(gel_markers_id, needle_id); % gel markers will correspond to next n_gel_markers largetst groups of voxels
totList = linspace(1, CC.NumObjects, CC.NumObjects);
stage_markers_id = setdiff(totList, [needle_id, gel_markers_id]); % stage markers will be the last few groups; it is set last because some markers might be invisible

%% Segment markers
% Voxel to MM scale
V2MM_scale = delMM*eye(3);

% Stage markers
stageMarkerVol = segmentVolume(size(newVol), CC, stage_markers_id);
stageMarkerCentroidsTab = regionprops3(stageMarkerVol, 'Centroid');
stage_marker_CT_coord= table2array(stageMarkerCentroidsTab)*V2MM_scale;
% volumeViewer(stageMarkerVol); % Visualize the volume


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

%% Registration

% Gel -> Find (R_gel, t_gel) to bring volume to design
[R_gel, t_gel, err_gel] = icp_wrapper(gel_marker_design_coord, gel_marker_CT_coord, 1);
gel_marker_CT_coord = apply_transform(gel_marker_CT_coord, R_gel, t_gel);

% Stage -> 1. Use (R_gel, t_gel) to transform the CT volume
stage_marker_CT_coord = apply_transform(stage_marker_CT_coord, R_gel, t_gel);

% Stage -> 2. Find (R_stage, t_stage) and bring design to volume
[R_stage, t_stage, err_stage] = icp_wrapper(stage_marker_CT_coord, stage_marker_design_coord, 1);
stage_marker_design_coord = apply_transform(stage_marker_design_coord, R_stage, t_stage);

%% Segment, Transform, and Estimate Needle Shape
needleMM = apply_transform(needleMM, R_gel, t_gel);
[y_fit_obj, y_gof] = fit(needleMM(:, 1), needleMM(:, 2), 'linearinterp');
[z_fit_obj, z_gof] = fit(needleMM(:, 1), needleMM(:, 3), 'linearinterp');
xp = (t_stage(1):1:max(needleMM(:, 1)))'; % wrt moved origin
yp = feval(y_fit_obj, xp);
zp = feval(z_fit_obj, xp);

Nel = length(xp);
% Nel = floor(length(xp)/2);

%% Potential error here...
% t_stage = [xp(1); yp(1); zp(1)];

%% Error Checking and Reporting
fprintf('Gel markers avg. registration error: %.3f mm\n', err_gel);
fprintf('Stage markers avg. registration error: %.3f mm\n', err_stage);

% Sanity checks
fprintf('Slice spacing: %.4f\n', delMM);
fprintf('Number of elements: %d\n', Nel);
fprintf('Size of each element: %.2f mm\n', (max(xp) - min(xp))/Nel);
fprintf('Origin x-axis separation: %.3f mm\n', t_stage(1));% x value should be close to the distance between two ames in CAD
fprintf('Insertion depth: %.2f mm\n', xp(end));
theta = atan2(R_stage(2, 1), R_stage(1, 1)); % using ZYX Euler angle, extracting only rotation about Z
fprintf('Baseline rotation angle (~0): %.3f degrees\n', theta/pi*180);% this is the rotation angle

%% Visualization
plotScene(gel_marker_design_coord, gel_marker_CT_coord, ...
   stage_marker_design_coord, stage_marker_CT_coord, ...
   [xp, yp, zp], needleMM, R_stage, t_stage);
title('Baseline Image')
end
