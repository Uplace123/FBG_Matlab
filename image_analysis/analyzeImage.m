function [xp, yp, zp] = analyzeImage(dataset, threshold)
% A heavily simplified version, disregarding stage data
%% Read CT
fprintf("Reading %s\n", dataset);
small_vol_threshold = 40;
[newVol, delMM] = readCT(dataset, threshold, small_vol_threshold);
% volumeViewer(newVol);

%% Read Body Frames
% Gel Frame
gel_marker_design_coord = readMarkerDesignCoord("gel_marker_coord.txt");
n_gel_markers = size(gel_marker_design_coord, 1);

%% Assign ID based on connected voxels
CC = bwconncomp(newVol); % find connected regions
numPixels = cellfun(@numel,CC.PixelIdxList);
[~, needle_id] = maxk(numPixels, 1); % needle will correspond to the largest group of connected voxels
[~, gel_markers_id] = maxk(numPixels, 1 + n_gel_markers);
gel_markers_id = setdiff(gel_markers_id, needle_id); % gel markers will correspond to next n_gel_markers largetst groups of voxels

%% Segment markers
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

%% Registration

% Gel -> Find (R_gel, t_gel) to bring volume to design
[R_gel, t_gel, err_gel] = icp_wrapper(gel_marker_design_coord, gel_marker_CT_coord, 1);
gel_marker_CT_coord = apply_transform(gel_marker_CT_coord, R_gel, t_gel);

%% Segment, Transform, and Estimate Needle Shape
needleMM = apply_transform(needleMM, R_gel, t_gel);
[y_fit_obj, y_gof] = fit(needleMM(:, 1), needleMM(:, 2), 'linearinterp');
[z_fit_obj, z_gof] = fit(needleMM(:, 1), needleMM(:, 3), 'linearinterp');
xp = (min(needleMM(:, 1)):.5:max(needleMM(:, 1)))'; % wrt moved origin
yp = feval(y_fit_obj, xp);
zp = feval(z_fit_obj, xp);
end
