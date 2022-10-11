addpath('./Galil_MATLAB_API/');
address = '192.168.1.201';

g = galil_connect(address); % A connection needs to be made before this line
res = galil_command(g, [char(18), char(22)]);
% disp(res)

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


flag = true;
promptX = "Input next X displacement (mm) negative-forward:";
promptY = "Input next Y displacement (mm) negative-left:";
while flag
    usr_input_X = input(promptX, 's');
    input_num_X = str2num(usr_input_X);
    
    usr_input_Y = input(promptY, 's');
    input_num_Y = str2num(usr_input_Y);
    tic
    if isempty(input_num_X) || isempty(input_num_Y)
        flag = false;
        disp('Quitting')
        continue
    else
        RelPosX = input_num_X;
        RelPosY = input_num_Y;
        Input_RelPos_X = round(RelPosX*1000);
        Input_RelPos_Y = round(RelPosY*1000);
        give_pos=strcat('PR ,,',num2str(Input_RelPos_Y),',',  num2str(Input_RelPos_X), ',');
        galil_command(g, give_pos);
        galil_command(g, 'BG CD');
        galil_command(g, 'MC CD');
    end
    toc
end


%% Homing before quitting
galil_command(g, 'PA ,,0,0,');
galil_command(g, 'BG CD');

%% Disconnect
galil_disconnect(g)
clear