clear;
clc;

% untested code
% run ini_interrogator.m
interrogator = ini_interrogator('IPaddress','192.168.1.11','Port',1852,'ReadTimeout',0.1);

RawData = [];

% get raw data
while 1
    RawData = Read_interrogator(1,2,4,interrogator);
    disp(RawData);
end

close_pnet();
