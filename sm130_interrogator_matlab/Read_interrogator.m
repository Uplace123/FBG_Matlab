function RawData = Read_interrogator(ReadCount,ChannelNumber,AANumber, varargin)

p =inputParser;
% optional parameters
addParameter(p,'IPaddress','192.168.1.11')
addParameter(p,'Port',1852)
addParameter(p,'ReadTimeout',0.1)
parse(p,varargin{:});

IPaddress = p.Results.IPaddress;
Port = p.Results.Port;
ReadTimeout = p.Results.ReadTimeout;

% ReadCount, quit the function after this number of times collect from
% interrogator.
% ChannelNumber, number of FBG
% AANumber, number of activate areas in one FBG.

% output, matrix
RawData = zeros(ReadCount,AANumber*ChannelNumber);

Interrogator = pnet('tcpconnect',IPaddress,Port);
% require tcp/udp/ip toolbox
% pnet('tcpconnect','hostname',port)
if (Interrogator < 0)
    disp("connect to interrogator failed!");
    return;
end

pnet(Interrogator,'setreadtimeout',ReadTimeout);
% con. unit second, blocks in 0.1s before timeout

pnet(Interrogator,'printf','#SET_CH_GAIN_DB 1 6.0 \n');
% print a formated string to the connection con

bytes_length = pnet(Interrogator,'Read',10,'char');
% pnet(con,'read',[size],[datatype],[swapping],'view','noblock')
% read an array of elements from a connection
% 'size' default 65536, return row vector if input is scalar.
% 'datatype' default 'char'
% pnet(con,'readline',[limitsize])
% read a string of characters until the newline character are reached.
% limitsiz, default 65536.
% reply from interogator about how many byte of total data 
% (normal should be 124 for 3 channels)

header = pnet(Interrogator,'Read',str2double(bytes_length),'char');
% read header

count = 1;

% tic
while 1

    pnet(Interrogator,'printf','#GET_UNBUFFERED_DATA \n');
    % command to get new set of data from the interogator

    bytes_length = pnet(Interrogator,'Read',10,'char');
    
    header = pnet(Interrogator,'Read',44,'int16','intel');
    %read header

    length_of_signal = (str2double(bytes_length)-88)/4;
    %get the length of signal

    current_data = zeros(1,AANumber*ChannelNumber);
  
    current_data = double(pnet(Interrogator,'Read',length_of_signal,'int32','intel'))/1e6;
    
    % try to get additional data
    test = double(pnet(Interrogator,'Read',1,'int32','intel'))/1e6;
    if ~isempty(test)
        disp("signal didn't fully read!");
        break;
    end
    
    RawData(count,:) = current_data;
    
    if count == ReadCount
        break;
    end

    count = count + 1;

end

pnet('closeall');
% close all pnet connections/sockets used in this matlab session

end