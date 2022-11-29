clear;
clc;

addpath ../sm130_interrogator_matlab/
% tested code
% run ini_interrogator.m
interrogator = ini_interrogator('IPaddress','192.168.1.11','Port',1852,'ReadTimeout',0.1);

RawData = [];

% get raw data
while 1
    %tic
    RawData = Read_interrogator(1,3,4,interrogator);
    disp(RawData);
    %toc % about 0.01s
end

close_pnet();
