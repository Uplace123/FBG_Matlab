clear;
clc;

% test

temp = zeros(1,8);
RawData = [];

% get 20 set of data
for i = [1:20]

    res_data = Read_interrogator(1,2,4,'IPaddress','192.168.1.11','Port',1852,'ReadTimeout',0.1);
    % RawData = Read_interrogator(ReadCount,ChannelNumber,AANumber,...)
    % RawData matrix(ReadCount,ChannelNumber * AANumber)
    RawData = [RawData; res_data];

end
disp(RawData);