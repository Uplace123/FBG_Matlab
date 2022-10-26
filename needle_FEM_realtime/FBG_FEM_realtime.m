function [ds, ks, xs] = FBG_FEM_realtime(sb, l, Db, Kb, ti, Nel, Mu, Alpha,...
                                        Interval,NumChannel,NumAA,interrogator, ...
                                        RefData,AA_lcn,FBG_switch)
% FBG_switch, 1 turn on FBG sensing, 0 turn off
if FBG_switch == 1
    % read rawdata
    RawData = Read_interrogator(1,NumChannel,NumAA,interrogator);
    % calculate curvautures base on rawdata
    curvatures = data_process(RawData,RefData,NumChannel,NumAA);
    curvatures_xy = curvatures(:,1);
    curvatures_xz = curvatures(:,2);

    % FEM outputs
    [ds, ks, xs] = FBG_integrated_FEM_Ogden_UTru(sb, l, Db, Kb, ti, Nel, Mu, ...
        Alpha, Interval, curvatures_xz, AA_lcn);
else
    % turn off FBG, data process
    curvatures_xz = [];
    [ds, ks, xs] = FBG_integrated_FEM_Ogden_UTru(sb, l, Db, Kb, ti, Nel, Mu, ...
        Alpha, Interval, curvatures_xz, AA_lcn);
end

end