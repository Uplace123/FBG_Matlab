addpath ./Galil_MATLAB_API/;
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

galil_command(g, 'KP ,,378,220,10');
galil_command(g, 'KI ,,145,84,6');
galil_command(g, 'KD ,,2831,1643,64');
galil_command(g, 'AC ,,5000,5000,50000');
galil_command(g, 'DC ,,5000,5000,50000');
galil_command(g, 'SP ,,8000,8000,100000');


flag = true;
promptX = "Input next X displacement (mm) negative-forward:";
promptY = "Input next Y displacement (mm) negative-left:";
promptK = "Input next k(Deg):";
while flag
    usr_input_X = input(promptX, 's');
    input_num_X = str2num(usr_input_X);
    
    usr_input_Y = input(promptY, 's');
    input_num_Y = str2num(usr_input_Y);
    
    usr_input_K = input(promptK, 's');
    input_num_K = str2num(usr_input_K);

    tic
    if isempty(input_num_X) || isempty(input_num_Y) || isempty(input_num_K)
        flag = false;
        disp('Quitting')
        continue
    else
        RelPosX = input_num_X;
        RelPosY = input_num_Y;
        RelPosK = input_num_K;
        Input_RelPos_X = round(RelPosX*1000);
        Input_RelPos_Y = round(RelPosY*1000);
        Input_RelPos_K = round(RelPosK*7031.25);
        give_pos=strcat('PR ,,',num2str(Input_RelPos_Y),',',  num2str(Input_RelPos_X), ',',num2str(Input_RelPos_K));
        galil_command(g, give_pos);
        galil_command(g, 'BG CDE');
        galil_command(g, 'MC CDE');
    end
    toc
end


%% Homing before quitting
galil_command(g, 'PA ,,0,0,0');
galil_command(g, 'BG CDE');

%% Disconnect
galil_disconnect(g)
clear