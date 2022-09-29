function Interrogator = ini_interrogator(varargin)
% initialize interrogator
% return tcp object

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

