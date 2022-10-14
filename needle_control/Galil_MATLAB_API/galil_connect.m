function [ gCon ] = galil_connect( ip_address )
%GALIL_CONNECT Summary of this function goes here
%   Detailed explanation goes here
    %gCon = tcpip(ip_address, 23, 'NetworkRole', 'client', 'Terminator', 'CR', 'Timeout', 12);
    gCon = tcpclient(ip_address, 23, 'Timeout', 12);
    fopen(gCon);
end

