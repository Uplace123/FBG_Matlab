function [ds, ks] = ogden_constrained_true(sb, l, Db, Kb, ti, Nel , Mu, Alpha)
% Flexible beam bending with Ogden1 hyperelastic material model with
% constrained imcompressibility condition. Solution using nonlinear FEM
% with load stepping
% Yanzhou Wang
% Dec 8 2021
addpath ../needle_FEM_realtime/WYZ_FEM/ 
%% Initialization
% Add Cholesky inverse
% addpath(genpath('./invChol/'));

% Constants
L = sb + l;
% Kb = -Db/(sb + rcm);


%FEM_params;
E = 250*1000;
OD = 1.27;
I = pi/4*(OD/2)^4; % in mm^4
ifplot = 1;


Interval = {[-sb, 0]; [0, l]};
MuT = [0; Mu]*10^-6; % Pa, but in mm^2; 1Pa = 1e-6 N/mm^2; 1kPa = 1e-3 N/mm^2
AlphaT = [0; Alpha]; % need abs(alpha) > 1 
GammaT = [0; 0];
PropertyTable = table(Interval, MuT, AlphaT, GammaT);
% ti = 25; % Initial length of undeformed tissue. Will affect solution

% FEM-specific constants
% Nel = 10; % Total umber of elements
Nen = 4; % Number of element DOF
DOF = 1:(2*Nel + 2);
nDOF = length(DOF);
d = zeros(nDOF, 1); % Initial guess of solution
h = L/Nel; % Finite element size
EBC = [1; 2]; % Displacement and slope of first element left node is prescribed
freeDOF = DOF;
freeDOF(EBC) = [];

% Algorithm-specific constants
max_inner_iter = 5; % Max number of iterations for Newton's method
max_outer_iter = 50; % Max number of iterations for load stepping
inner_iter = 0; % Initialize inner_iter counter
outer_iter = 0; % Initialize outer_iter counter
tol = 1e-3; % Convergence criterion
converged = 0;
load_ratio = 1;
EBC_cur = zeros(2, 1); % For load stepping
EBC_converged = zeros(2, 1); % For load stepping from previously converged EBC

% Construction of LM
LM = zeros(Nen, Nel);
for e = 1:Nel
    LM(1, e) = 2*e - 1;
    LM(2, e) = 2*e - 0;
    LM(3, e) = 2*e + 1;
    LM(4, e) = 2*e + 2;
end

%% FEM Main
% tic
% Load stepping
while (~all(EBC_cur == EBC)) && (outer_iter < max_outer_iter)
    outer_iter = outer_iter + 1;
    inner_iter = 0;
    converged = 0;
    EBC_cur = EBC_converged + load_ratio*EBC;
    % FEM solution
    while inner_iter < max_inner_iter
        K = zeros(nDOF, nDOF);
        F = zeros(nDOF, 1);
        P = zeros(nDOF, 1);
        % Apply EBC
        d(EBC) = [Db; Kb];
        % Loop over each element
        for e = 1:Nel
            % Get local nodal values from global d
            d_i_local = d(2*e - 1:2*e + 2);
            [ke, pe] = compute_element_matrix(d_i_local, ti, E, I, PropertyTable, e, h, -sb);
            % Global assembly process
            for i = 1:4
                for j = 1:4
                    K(LM(i, e), LM(j, e))=K(LM(i, e), LM(j, e)) + ke(i, j);
                end
                P(LM(i, e))=P(LM(i, e))+pe(i);
            end
        end

        % Newton's method
        % Use only free DOF from the list of DOF to compute d
        delta_d = pinv(K(freeDOF, freeDOF))*(F(freeDOF, 1)-P(freeDOF, 1));
        
        if(norm((F(freeDOF, 1)-P(freeDOF, 1))) <= tol)
            converged = 1;
            %fprintf("Converged at inner step %d, outer step %d\n", inner_iter + 1, outer_iter);
            break;
        end
        d(freeDOF, 1) = d(freeDOF, 1) + delta_d;
        inner_iter = inner_iter + 1;
    end

    % Check for convergence of current inner iteration
    if converged ~= 1
        load_ratio = 0.5*load_ratio;
        EBC_cur = zeros(2, 1);
        % fprintf("No convergence. Decreasing load step\n");
    else
        EBC_converged = EBC_cur;
    end
end

% output only the displacement
ds = extract_dist_from_d(d);
ks = extract_slop_from_d(d);

%% Display of final result
% if converged
%     %disp('Newton-Ralphson converged');
% else
%     %disp('Newton-Ralphson does not converge')
% end
% % toc
% 
% % Plot result
% if ifplot == 1 && converged
%     plot_result(ds, sb, h, l, Nel, PropertyTable);
% end
% 
% end

%% Defined helper functions
% Look up material properties at given location
function [MuT_e, AlphaT_e, GammaT_e] = lookup_property(PropertyTable, cur_e, ele_size, x_begin)
x_mid = x_begin + cur_e*ele_size - 1/2*ele_size; % Find the midpoint location of the element in global coordinate
for i = 1:numel(PropertyTable.Interval)
    cur_interval = PropertyTable.Interval{i};
    if (x_mid >= cur_interval(1)) && (x_mid < cur_interval(2))
        MuT_e = PropertyTable.MuT(i);
        AlphaT_e = PropertyTable.AlphaT(i);
        GammaT_e = PropertyTable.GammaT(i);
    end
end
end

% Compute element stiffness matrix and internal force vector
function [ke, pe] = compute_element_matrix(d_i_local, ti, E, I, PropertyTable, e, h, x_begin)
% Find properties
[MuT_e, AlphaT_e, GammaT_e] = lookup_property(PropertyTable, e, h, x_begin);
% Calculate integrals
% _beam will stay the same
pe_beam = calc_pe_beam(d_i_local, h, E, I);
ke_beam = calc_ke_beam(h, E, I);
% _cont will depend on external force formulation
pe_cont = calc_pe_cont(d_i_local, h, ti, MuT_e, AlphaT_e, GammaT_e);
ke_cont = calc_ke_cont(d_i_local, h, ti, MuT_e, AlphaT_e, GammaT_e);

pe = pe_beam + pe_cont;
ke = ke_beam + ke_cont;
end


%% Integration subroutines
function pe_beam = calc_pe_beam(d_i_local, h, E, I)
% Element internal force vector of the beam
pe_beam = calc_ke_beam(h, E, I)*d_i_local;
end

function ke_beam = calc_ke_beam(h, E, I)
% Element stiffness from beam bending
ke_beam = E*I/h^3*[12,   6*h,    -12,    6*h;
    6*h,  4*h^2,  -6*h,   2*h^2;
    -12,  -6*h,   12,     -6*h;
    6*h,  2*h^2,  -6*h,   4*h^2];
end


function pe_cont = calc_pe_cont(d_i_local, h, ti, MuT_e, AlphaT_e, GammaT_e)
% Element internal force vector from tissue contact
ig = [-1/sqrt(3); 1/sqrt(3)]; % Gauss integration points
wg = [1; 1]; % Gauss weights
N = [shape(ig(1), h); shape(ig(2), h)]; % Shape function at each gauss point
N_zeta = [dshape(ig(1), h); dshape(ig(2), h)]; % Shape function derivative at each gauss point
u = [N(1, :)*d_i_local; N(2, :)*d_i_local]; % Displacements at each gauss points
u_zeta = [N_zeta(1, :)*d_i_local; N_zeta(2, :)*d_i_local];% Displacements derivative at each gauss points

% Sum over all gauss points
pe_cont = zeros(1, 4);
for i = 1:length(ig)
    pe_cont = pe_cont + ...
        wg(i)*(...
        N(i, :)*2*MuT_e*(((ti - abs(u(i)))/ti)^(AlphaT_e - 1) + ((ti - abs(u(i)))/ti)^(-AlphaT_e - 1))*...
        u(i)*(1 - GammaT_e*sin(atan(u_zeta(i)*(2/h)))^2)*(h/2)...
        );
end
pe_cont = pe_cont'; % transpose to get the right dimension
end

function ke_cont = calc_ke_cont(d_i_local, h, ti, MuT_e, AlphaT_e, GammaT_e)
% Element stiffness from tissue contact
ig = [-1/sqrt(3); 1/sqrt(3)]; % Gauss integration points
wg = [1; 1]; % Gauss weights
N = [shape(ig(1), h); shape(ig(2), h)]; % Shape function at each gauss point
N_zeta = [dshape(ig(1), h); dshape(ig(2), h)]; % Shape function derivative at each gauss point
u = [N(1, :)*d_i_local; N(2, :)*d_i_local]; % Displacements at each gauss points
u_zeta = [N_zeta(1, :)*d_i_local; N_zeta(2, :)*d_i_local];% Displacements derivative at each gauss points

% Sum over all gauss points
ke_cont = zeros(4, 4);
for i = 1:length(ig)
    ke_cont = ke_cont + ...
        wg(i)*(...
        N(i, :)'*2*MuT_e*((AlphaT_e - 1)*((ti - abs(u(i)))/ti)^(AlphaT_e - 2)*(-1/ti*sign(u(i))*N(i, :)) + ...
        (-AlphaT_e - 1)*((ti - abs(u(i)))/ti)^(-AlphaT_e - 2)*(-1/ti*sign(u(i))*N(i, :)))*u(i)* ...
        (1 - GammaT_e*sin(atan(u_zeta(i)*(2/h)))^2)*(h/2) + ...
        ...
        N(i, :)'*2*MuT_e*(((ti - abs(u(i)))/ti)^(AlphaT_e - 1) + ((ti - abs(u(i)))/ti)^(-AlphaT_e - 1))*N(i, :)*...
        (1 - GammaT_e*sin(atan(u_zeta(i)*(2/h)))^2)*(h/2) + ...
        ...
        N(i, :)'*2*MuT_e*(((ti - abs(u(i)))/ti)^(AlphaT_e - 1) + ((ti - abs(u(i)))/ti)^(-AlphaT_e - 1))*u(i)*...
        -GammaT_e*((2*u_zeta(i)*(2/h))/((u_zeta(i)*(2/h))^2 + 1)^2*N(i, :)*(2/h) - (2*(u_zeta(i)*(2/h))^3)/((u_zeta(i)*(2/h))^2 + 1)^2*N(i, :)*(2/h))*(h/2)...
        );
end
end

%% Shape functions
function N = shape(zeta, h)
% Shape functions N(zeta)
N = [(zeta -1)^2*(zeta + 2)/4, ...
    h*(zeta - 1)^2*(zeta + 1)/8, ...
    (zeta + 1)^2*(zeta + 2)/4, ...
    h*(zeta - 1)*(zeta + 1)^2/8];
end

function N_zeta = dshape(zeta, h)
% Shape functions N_zeta(zeta)
N_zeta = [3*(zeta - 1)*(zeta + 1)/4, ...
    h*(zeta - 1)*(3*zeta + 1)/8, ...
    3*(zeta + 1)*(1 - zeta)/4, ...
    h*(zeta + 1)*(3*zeta - 1)/8];
end

function N_zeta_zeta = ddshape(zeta, h)
% Shape functions N_zeta_zeta(zeta)
N_zeta_zeta = [3*zeta/2, ...
    (3*zeta - 1)*h/4, ...
    -3/2*zeta, ...
    (3*zeta + 1)*h/4];
end


% %% Defined auxiliary functions
% function plot_result(d, sb, h, l, Nel, PropertyTable)
% % Plotting needle shape
% x = -sb:h:l;
% figure('Name', sprintf('FEM Solution with Load Stepping'));
% plt = plot(x, rand(1, size(x, 2)), 'k-');
% % ds = zeros(Nel + 1, 1);
% % for i = 1:(Nel + 1)
% %     ds(i) = d(2*i - 1) ;
% % end
% grid on;
% set(plt, 'YData', d);
% axis equal;
% 
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
% end

end