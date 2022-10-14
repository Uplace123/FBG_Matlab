function [dx_c,dy_c,dr,Db,Kb] = robot_geometric(dy,dk,Db,Kb)
% this script get the dx dy at robot control position 
% base on dy and dk at FEM base position 
% this script also update Db in mm Kb in DEG 
% dy, displacement at y axis in mm
% dk, rotation in degree

D = 50.62; % in mm
alpha = atan(dk);
alpha_kb = atan(Kb);
x1 = D*cos(alpha_kb);
y1 = D*sin(alpha_kb);
% get current base location and slope
Db = Db + dy;
Kb = Kb + dk;
alpha_kb = atan(Kb);
%alpha_kb = alpha + alpha_kb;
x2 = D*cos(alpha_kb);
y2 = D*sin(alpha_kb);
dx_c = x1 - x2;
dy_c = y1 - y2 + dy;

dr = rad2deg(alpha);


end