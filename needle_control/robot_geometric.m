function [dx_c,dy_c,dk] = robot_geometric(dy,dk)
% this script get the dx dy at robot control position 
% base on dy and dk at FEM base position 
% dy, displacement at y axis
% dk, rotation

D = 50.62;
alpha = atan(deg2rad(dk));
dx_c = D * (1 - sin(pi/2 - alpha));
dy_c = dy + D * cos(pi/2 - alpha);

end