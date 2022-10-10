function segVol = segmentVolume(V_size, CC, idx)
segVol = false(V_size);
for i = 1:length(idx)
    segVol(CC.PixelIdxList{idx(i)}) = 1;
end
end