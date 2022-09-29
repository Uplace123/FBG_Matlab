% Plotting needle shape

% get size
[ds, ks, sb, h, l, PropertyTable] = FBG_integrated_FEM_Ogden_UTru([0;0;0;0]);
sz_ds = size(ds);
sz_ks = size(ks);
sz_sb = size(sb);
sz_h = size(h);
sz_l = size(l);

filename = fullfile(tempdir,'communicate.dat');

a = exist(filename);

if a == 0
    fid = fopen(filename,'w');
    fwrite(fid,ds,'double');
    fwrite(fid,ks,'double');
    fwrite(fid,sb,'double');
    fwrite(fid,h,'double');
    fwrite(fid,l,'double');
    fclose(fid);
end

m = memmapfile(filename, 'Writable',true, 'Format', ...
    {'double', sz_ds, 'ds';
    'double', sz_ks, 'ks';
    'double', sz_sb, 'sb';
    'double', sz_h, 'h';
    'double', sz_l, 'l';
    });


x = -m.Data.sb:m.Data.h:m.Data.l;
f = gcf;
set(f, 'Name', sprintf('FEM Solution with Load Stepping'));
plt = plot(x, rand(1, size(x, 2)), 'k-', 'Parent', gca);
grid on;
axis equal;
% plot some reference curve
% todo

% plot the curvature at AA
% todo

% plot the position at tip
% todo

while 1
    set(plt, 'YData', m.Data.ds);
    %disp(m.data.ds);
    drawnow
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
