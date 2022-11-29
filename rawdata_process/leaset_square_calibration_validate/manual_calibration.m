% pnetPath = pnetPath(:,1:end-20);
addpath ./sm130_interrogator_matlab/


% 6. Define the number of datapoints needed for each recording location
FBG_data_point = 200;
% 7. Define the number of repeatitions during validation movement
val_repetitions = 5;
% 8. Define number of active areas
channels = 2;
AAs = 4;
% 9. Define number of slots on calibration jig
cal_curve = [0,0.25,0.8,1.0,1.25]';


slot_num = length(cal_curve);

cal_rot = [0,90];
cal_avg = zeros(length(cal_rot),6,channels*AAs);
%path = strcat(cd,'/data/093022');
%mkdir (path);
calPath = cd;
mkdir(calPath);
for k = 1:5
    for i = 1:slot_num
        for n = 1:length(cal_rot)
        
            txt = strcat('This is slot no.',num2str(i),', ',num2str(cal_rot(n)), ' deg, trial ',num2str(k));
            disp(txt);
            disp('Press Enter when ready to collect data');
            pause;
            
            cal_data = fbgdata_collection(FBG_data_point,2,4);
            fileName = strcat(calPath,'/validation.xls');
            sheetName = strcat('trial',num2str(k),'_',num2str(cal_curve(i)),'mm','_',num2str(cal_rot(n)),'deg');
            writematrix(cal_data,fileName,'sheet',sheetName);
            
            cal_avg(n,i,:) = squeeze(cal_avg(n,i,:))' + mean(cal_data,1);
        end
        
    end
    
    cal_avg = cal_avg./k;
end


% cal_AA = {cal_avg(:,:,[1,3,5,7]),cal_avg(:,:,[2,4,5,8])};
% 
% for i = 1:AAs
%     % figure for calibration
%     figure(4)
%     subplot(1,2,1)
%     hold on
%     for j = 1:channels
%         plot([cal_curve;-cal_curve]',[cal_AA{1,i}(1,:,j),cal_AA{1,i}(3,:,j)],'*');
%     end
%     hold off
%     xlabel('curvature x (1/m)');
%     ylabel('T corr. signal response (nm)');
%     legend({'CH1','CH2','CH3'},'Location','eastoutside')
%     
%     subplot(1,2,2)
%     hold on
%     for j = 1:channels
%         plot([cal_curve;-cal_curve],[cal_AA{1,i}(2,:,j),cal_AA{1,i}(4,:,j)],'*');
%     end
%     hold off
%     xlabel('curvature y (1/m)');
%     ylabel('T corr. signal response (nm)');
%     legend({'CH1','CH2','CH3'},'Location','eastoutside')
%     
%     title(strcat('AA',num2str(i),' response vs curvature'),'left: x  right: y');
%     picName = strcat(calPath,'\calibration_AA',num2str(i),'.png');
%     saveas(gcf,picName);
%     clf(gcf);
%     
% end