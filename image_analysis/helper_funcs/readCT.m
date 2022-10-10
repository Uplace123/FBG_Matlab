function [newVol, delMM] = readCT(dataset, threshold, varargin)
% path = fullfile('~/Desktop', 'DICOM', dataset + '/'); % moved DICOM folder to desktop due to insufficient space in MATLAB Drive
[V, spatial] = dicomreadVolume(dataset);
V = squeeze(V);
newVol = V > threshold;
newVol = bwmorph3(newVol, 'majority');
if isempty(varargin)
    small_vol_threshold = 0;
else
    small_vol_threshold = varargin{1};
end
newVol = bwareaopen(newVol, small_vol_threshold);
delMM = unique(spatial.PixelSpacings);

if numel(delMM) ~= 1
    s1 = sprintf('Non-uniform slice spacing:');
    s2 = sprintf(repmat('%.4f ', 1, numel(delMM)), delMM);
    s3 = sprintf('\nUsing %.4f\n', delMM(1));
    warning(strcat(s1, s2, s3));
    delMM = delMM(1);
end

% flip the volume LR if necessary (check volume with actual setup)
% tformLR = affine3d([-1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1]);
% newVol = imwarp(newVol, tformLR);
end
