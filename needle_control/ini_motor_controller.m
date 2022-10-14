function g = ini_motor_controller(address)
g = galil_connect(address); % A connection needs to be made before this line
res = galil_command(g, [char(18), char(22)]);
% parameter of controller
% axis B: rotery stage                7031.25 unit: 1 degree
% axis C: Y direction of XYZ robot    1000 unit: 1 mm
% axis D: Z direction of XYZ robot    1000 unit: 1 mm
% axis E: X direction of XYZ robot    1000 unit: 1 mm
% galil_command(g, 'KP,10,54,54,54');
% galil_command(g, 'KI,6,4,4,4');
% galil_command(g, 'KD,64,480,480,480');
% galil_command(g, 'AC,5000,3000,3000,3000');
% galil_command(g, 'DC,5000,3000,3000,3000');
% galil_command(g, 'SP,30000,5000,5000,5000');

galil_command(g, 'KP ,,378,220,10');
galil_command(g, 'KI ,,145,84,6');
galil_command(g, 'KD ,,2831,1643,64');
galil_command(g, 'AC ,,5000,5000,50000');
galil_command(g, 'DC ,,5000,5000,50000');
galil_command(g, 'SP ,,8000,8000,100000');

end