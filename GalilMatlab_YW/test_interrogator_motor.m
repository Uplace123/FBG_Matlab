% test galil motor controller and interrogator at the same time

addpath ./Galil_MATLAB_API/
addpath ../sm130_interrogator_matlab/
% Galil controller ip address
address = '192.168.1.201';
g = galil_connect(address); % A connection needs to be made before this line
res = galil_command(g, [char(18), char(22)]);
%disp(res)
% tested code
% run ini_interrogator.m
interrogator = ini_interrogator('IPaddress','192.168.1.11','Port',1852,'ReadTimeout',0.1);

RawData = [];

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

galil_command(g, 'KP ,,378,220,');
galil_command(g, 'KI ,,145,84,');
galil_command(g, 'KD ,,2831,1643,');
galil_command(g, 'AC ,,3000,3000,');
galil_command(g, 'DC ,,3000,3000,');
galil_command(g, 'SP ,,5000,5000,');

RelPosX = 0.01;
RelPosY = 0;

for i = 1:1:50
    %tic
    tic
    Input_RelPos_X = round(RelPosX*1000);
    Input_RelPos_Y = round(RelPosY*1000);
    give_pos=strcat('PR ,,',num2str(Input_RelPos_Y),',',  num2str(Input_RelPos_X), ',');
    galil_command(g, give_pos);
    galil_command(g, 'BG CD');
    galil_command(g, 'MC CD');
    toc
    % get raw data
    RawData = Read_interrogator(1,2,4,interrogator);
    
    %disp(RawData);
end


%% Homing before quitting
galil_command(g, 'PA ,,0,0,');
galil_command(g, 'BG CD');

%% Disconnect
galil_disconnect(g)
close_pnet();
clear