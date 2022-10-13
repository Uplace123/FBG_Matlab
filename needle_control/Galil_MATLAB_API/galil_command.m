function [ result ] = galil_command( gCon, command )
%GALIL_COMMAND Summary of this function goes here
%   Detailed explanation goes here
    cmd = sprintf('%s\r', command);
    % disp(cmd);
    fwrite(gCon, cmd);
    while gCon.BytesAvailable==0
    end
    if  gCon.BytesAvailable~=0
        data = fread(gCon, gCon.BytesAvailable);
        %str = cast(data', 'char');
        
        response = strsplit(cast(data', 'char'), ':');
        result = response{1,1};
    else 
        %disp('Nothing to read');
    end
end

