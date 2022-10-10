function [ds, ks, xs] = FBG_FEM_realtime(sb, l, Db, Kb, ti, Nel, Mu, Alpha, Interval,NumChannel,NumAA,interrogator,RefData,AA_lcn)

% read rawdata
RawData = Read_interrogator(1,NumChannel,NumAA,interrogator);
% calculate curvautures base on rawdata
curvatures = data_process(RawData,RefData,NumChannel,NumAA);
curvatures_xy = curvatures(:,1);
curvatures_xz = curvatures(:,2);
% FEM outputs
[ds, ks, xs] = FBG_integrated_FEM_Ogden_UTru(sb, l, Db, Kb, ti, Nel, Mu, ...
    Alpha, Interval, curvatures_xz, AA_lcn);

end