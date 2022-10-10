function marker_coord_mat = readMarkerDesignCoord(filename)
%% Set File Read Options
opts = delimitedTextImportOptions("NumVariables", 3);
% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = ",";
% Specify column names and types
opts.VariableNames = ["VarName1", "VarName2", "VarName3"];
opts.VariableTypes = ["double", "double", "double"];
% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
%% Read Body Frames
marker_coord_tab = readtable(filename, opts);
marker_coord_mat = table2array(marker_coord_tab);
end