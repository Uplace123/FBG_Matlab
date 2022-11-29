function fbg_data = fbgdata_collection(FBG_data_point,num_CH,num_AA) % this is for 3CH 4AA
% FBG data data collection from interrogator and Shape Sensing
% tic

%fbg_data = zeros(FBG_data_point,9);
num_datapts = num_CH*num_AA;
fbg_data = zeros(FBG_data_point,num_datapts);

Interrogator_t = pnet('tcpconnect','192.168.1.11',1852);
pnet(Interrogator_t,'setreadtimeout',0.1)
pnet(Interrogator_t,'printf','#SET_CH_GAIN_DB 1 6.0 \n');

bytes_length = pnet(Interrogator_t,'Read',10,'char');%reply from interogator about how many byte of total data (normal should be 124 for 3 channels)
pnet(Interrogator_t,'Read',str2double(bytes_length),'char');%read header


% tic
for j=1:FBG_data_point
    %every 9 times loop start sending command to get new set of data from the interogator
    pnet(Interrogator_t,'printf','#GET_UNBUFFERED_DATA \n');%command to get new set of data from the interogator
    %pnet(Interrogator_t,'Readline');%read header
    bytes_length = pnet(Interrogator_t,'Read',10,'char');%reply from interogator about how many byte of total data (normal should be 124 for 3 channels)
    pnet(Interrogator_t,'Read',44,'int16','intel');%read header
    length_of_signal = (str2double(bytes_length)-88)/4;
    current_data = zeros(1,num_datapts+1);
    %disp(length_of_signal)
    for i = 1:length_of_signal
        current_data(1,i) = double(pnet(Interrogator_t,'Read',1,'int32','intel'))/1e6;
    end
    
    if length_of_signal == 10
        current_data(6) = [];
    else
        current_data(end) = [];
    end
    
    fbg_data(j,:) = current_data;
end

pnet('closeall');

end