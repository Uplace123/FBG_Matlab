function RawData = Read_interrogator(ReadCount,ChannelNumber,AANumber, Interrogator)

% output, matrix
RawData = zeros(ReadCount,AANumber*ChannelNumber);

count = 1;
% tic
for i = 1:ReadCount

    pnet(Interrogator,'printf','#GET_UNBUFFERED_DATA \n');
    % command to get new set of data from the interogator

    bytes_length = pnet(Interrogator,'Read',10,'char');
    
    header = pnet(Interrogator,'Read',44,'int16','intel');
    %read header
    length_of_signal = (str2double(bytes_length)-88)/4;
    %get the length of signal

    current_data = double(pnet(Interrogator,'Read',length_of_signal,'int32','intel'))/1e6;
    
    % try to get additional data
    %test = double(pnet(Interrogator,'Read',1,'int32','intel'))/1e6;
    
    
%     if ~isempty(test)
%         disp("signal didn't fully read!");
%         break;
%     end
    
    RawData(count,:) = current_data;
    
    count = count + 1;
end


end