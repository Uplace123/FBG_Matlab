% Plotting needle shape
close all;
clear;
clc;

addpath ../needle_FEM_realtime/WYZ_FEM/
addpath ../needle_FEM_realtime/WYZ_FEM/invChol/

%% FEM needed input
load plot_params.mat

NumAA_real = 4; % keep the same as needle
AA_lcn = AA_lcn_base;
NumAA = size(AA_lcn,2); % this is the number off AA used in FEM

% get size
curvatures = zeros(NumAA_real,2);
curvatures_xz = curvatures(:,2);
% position of AA, dist from the baseS
tip_lcn = l + sb;


[ds, ks, xs] = FBG_integrated_FEM_Ogden_UTru(sb, l, Db , Kb, ti, Nel, Mu, ...
                                             Alpha, Interval, curvatures_xz, AA_lcn);

sz_ds = size(ds);
sz_ks = size(ks);
sz_xs = size(xs);
sz_curvature = size(curvatures);

filename = fullfile(tempdir, 'communicate.dat');

a = exist(filename);

if a == 0
    fid = fopen(filename, 'w');
    fwrite(fid, ds, 'double');
    fwrite(fid, ks, 'double');
    fwrite(fid, xs, 'double');
    fwrite(fid, curvatures, 'double');
    fclose(fid);
end

m = memmapfile(filename, 'Writable',true, 'Format', ...
    {'double', sz_ds, 'ds';
    'double', sz_ks, 'ks';
    'double', sz_xs, 'xs';
    'double', sz_curvature, 'curvature';
    });



% x = -m.Data.sb:m.Data.h:m.Data.l;

% get the ydata at AA

index_aa  = (AA_lcn  / (l+sb))*(size(xs,2) - 1);
index_tip = (tip_lcn / (l+sb))*(size(xs,2) - 1);


f = gcf;
set(f, 'Name', sprintf('FEM Solution with Load Stepping'));
set(f, 'position', [200,200,2000,1500]);
plt1 = plot(xs, rand(1, size(xs, 2)), 'k-', 'Parent', gca,'LineWidth',5);
hold on
plt2 = plot(xs(index_aa+1), rand(1, size(index_aa, 2)), 'r*', 'Parent', gca,'LineWidth', 5);
hold on
plt3 = plot(xs(index_tip+1), rand(1, size(index_tip, 2)), 'b*', 'Parent', gca,'LineWidth', 5);
hold on

% plot reference curves
% offset = 35;
% load refer_curve_05.mat
% plot(x_1*1000 + offset,y_1*1000,'--');
% load refer_curve_16.mat
% plot(x_1*1000 + offset,y_1*1000,'--');
% load refer_curve_20.mat
% plot(x_1*1000 + offset,y_1*1000,'--');
% load refer_curve_25.mat
% plot(x_1*1000 + offset,y_1*1000,'--');
% load refer_curve_32.mat
% plot(x_1*1000 + offset,y_1*1000,'--');


tip_text = text(0,0,"",'FontSize',25,'Parent',gca);
AA1_text = text(10,-40,"",'FontSize',25,'Parent',gca);
AA2_text = text(40,-40,"",'FontSize',25,'Parent',gca);
AA3_text = text(70,-40,"",'FontSize',25,'Parent',gca);
AA4_text = text(100,-40,"",'FontSize',25,'Parent',gca);


grid on

xlabel("length",'FontSize',30);
ylabel("deflection",'FontSize',30);
set(gca, "FontSize", 30);

axis equal;
xlim([-sb, l+20]);
ylim([-60,60]);
% plot some reference curve
% todo


while 1
    %tic
    set(plt1, 'YData', m.Data.ds);
    
    Y_AA = m.Data.ds(index_aa + 1);
    Y_tip = m.Data.ds(index_tip + 1);
    % plot AA area and tip
    set(plt2, 'YData', Y_AA);
    set(plt3, 'YData', Y_tip);


    % tip deflection
    str = "Tip deflection: "+ newline + num2str(Y_tip)+"mm";
    %tip_text = text(index_tip-2,Y_tip+2,cellstr(str),'FontSize',15);
    set(tip_text,'Position',[xs(index_tip+1)+2,Y_tip], 'String', str);
    % curvature at AA
    for i = 1:NumAA
        switch i
            case 1 
                str_A1 = "Curvature AA1: " + newline +num2str(m.Data.curvature(1,2))+" 1/m";
                set(AA1_text,'Position',[xs(index_aa(1)+1)-10,Y_AA(1)+10], 'String', str_A1);
            case 2
                str_A2 = "Curvature AA2: " + newline +num2str(m.Data.curvature(2,2))+" 1/m";
                set(AA2_text,'Position',[xs(index_aa(2)+1)-10,Y_AA(2)-10], 'String', str_A2);
            case 3
                str_A3 = "Curvature AA3: " + newline +num2str(m.Data.curvature(3,2))+" 1/m";
                set(AA3_text,'Position',[xs(index_aa(3)+1)-10,Y_AA(3)+10], 'String', str_A3);
            case 4
                str_A4 = "Curvature AA4: " + newline +num2str(m.Data.curvature(4,2))+" 1/m";
                set(AA4_text,'Position',[xs(index_aa(4)+1)-10,Y_AA(4)-10], 'String', str_A4);
        end
    end

    drawnow
    pause(0.05);
    %toc
end


% % Plotting patches
% intervals = PropertyTable.Interval;
% num_inter = size(intervals, 1);
% Ps = [];
% titles = [];
% Py_Lims = get(gca, 'YLim');
% for j = 1:num_inter
%     Px = [intervals{j}(1), intervals{j}(2), intervals{j}(2), intervals{j}(1)];
%     Py = [Py_Lims(1), Py_Lims(1), Py_Lims(2), Py_Lims(2)];
%     if j == 1
%         Ps = patch(Px, Py, 'w', 'FaceAlpha', 0.2);
%     else
%         Ps = [Ps, patch(Px, Py, j/num_inter, 'FaceAlpha', 0.2)];
%     end
%     titles = [titles; sprintf("Mu_{T%d} = %.0f Pa, Alpha_{T%d} = %.2f, Gamma_{T%d} = %.2f", ...
%         j - 1, PropertyTable.MuT(j)*10^6, ...
%         j - 1, PropertyTable.AlphaT(j), ...
%         j - 1, PropertyTable.GammaT(j))];
% end
% legend(Ps, titles);
% axis tight
