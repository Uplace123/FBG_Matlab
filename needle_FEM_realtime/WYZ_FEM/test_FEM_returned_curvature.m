%% FBG integrated FEM Ogden UTru test
% Provided that the FEM agrees with the ground truth CT scan, use converged
% FEM solution to reverse-calculate local curvature information at given AA
% locations
% The calculated curvatures are stored as a variable 'curvatures'

close all
%% Initialization
sb = 5;
l = 40;
Db = 5;
Kb = 0;
ti = 25;
Nel = sb + l;
Mu = 5e3;
Alpha = 5;
E = 250*1000; % 200GPa but in mm^2
OD = 1.27;
I = pi/4*(OD/2)^4; % in mm^4
ifplot = 1;

% Add Cholesky inverse
addpath(genpath('./invCol/'))

% Constants
L = sb + l;
Interval = {[-sb, 0]; [0, l]};
MuT = [0; Mu*10^-6]; % Pa, but in mm^2; 1Pa = 1e-6 N/mm^2; 1kPa = 1e-3 N/mm^2
AlphaT = [-2; Alpha]; % need abs(alpha) > 1
GammaT = zeros(size(MuT));
PropertyTable = table(Interval, MuT, AlphaT, GammaT);

Nen = 4; % Number of element DOF
DOF = 1:(2*Nel + 2);
nDOF = length(DOF);
d = zeros(nDOF, 1); % Initial guess of solution
h = L/Nel; % Finite element size
EBC_idx = [1; 2]; % Displacement and slope of first element left node is prescribed
freeDOF = DOF;
freeDOF(EBC_idx) = [];

% Algorithm-specific constants
max_inner_iter = 15; % Max number of iterations for Newton's method
max_outer_iter = 10; % Max number of iterations for load stepping
inner_iter = 0; % Initialize inner_iter counter
outer_iter = 0; % Initialize outer_iter counter
tol = 1e-5; % Convergence criterion
converged = 0;
load_ratio = 1;
EBC_des = [Db; Kb];
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

% FBG AA locations
AA_lcn = 1:2:Nel - 2;
AA_er = round(AA_lcn./h) + 1; % elements where the left-moment is fixed
AA_crv = AA_crv(AA_er >= 0); % curvatures that have negative element indices are skipped
AA_er = AA_er(AA_er >= 0); % elements that have negative indices are skipped

%% FEM Main
% Load stepping
while converged == 0 && (outer_iter < max_outer_iter)
    outer_iter = outer_iter + 1;
    inner_iter = 0;
    converged = 0;
    EBC_cur = EBC_converged + load_ratio*EBC_des;
    % FEM solution
    while inner_iter < max_inner_iter
        K = zeros(nDOF, nDOF);
        F = zeros(nDOF, 1);
        P = zeros(nDOF, 1);
        % Apply EBC
        d(EBC_idx) = EBC_cur;
        % Loop over each element
        for e = 1:Nel
            % Get local nodal values from global d
            d_i_local = d(2*e - 1:2*e + 2);
            [ke, pe] = compute_element_matrix(d_i_local, ti, E, I, ...
                PropertyTable, e, h, -sb);

            % Global assembly process
            % Using vectors instead of FOR loops
            K(LM(1:4, e), LM(1:4, e)) = K(LM(1:4, e), LM(1:4, e)) + ke;
            P(LM(1:4, e)) = P(LM(1:4, e)) + pe;
        end

        % Newton's method
        % Use only free DOF from the list of DOF to compute d
        delta_d = invChol_mex(K(freeDOF, freeDOF))*(F(freeDOF, 1)-P(freeDOF, 1));

        if(norm((F(freeDOF, 1)-P(freeDOF, 1))) <= tol && all(abs(EBC_cur - EBC_des) <= tol))
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
        disp("No convergence. Decreasing load step");
    else
        EBC_converged = EBC_cur;
        % disp("Converged. Continuing load stepping");
    end
end

% output only the displacement
ds = extract_dist_from_d(d);
ks = extract_slop_from_d(d);
dds = diff(ds);
xs = zeros(1, Nel + 1);
xs(1) = -sb;
for i = 2:(Nel + 1)
    xs(i) = xs(i - 1) + sqrt(h^2 - dds(i - 1)^2);
end

%% Display of final result
if converged
    disp('Newton-Ralphson converged');
else
    disp('Newton-Ralphson does not converge')
end

% Plot result
if ifplot == 1 && converged
    plot_result(ds, xs, PropertyTable);
end

%% Calculate curvatures at AA locations based on converged FEM
if converged
    curvatures = zeros(size(AA_er));
    for i = 1:length(AA_er)
        d_AA_local = d(2*AA_er(i) - 1:2*AA_er(i) + 2);
        pe_beam = calc_ke_beam(h, E, I)*d_AA_local;
        curvatures(i) = pe_beam(2)/E/I;
    end
end

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
function [ke, pe] = compute_element_matrix(d_i_local, ti, E, I, ...
    PropertyTable, e, h, x_begin)
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
        N(i, :)*2*MuT_e*(((ti - abs(u(i)))/ti)^(AlphaT_e - 1) + 1/2*((ti - abs(u(i)))/ti)^(-AlphaT_e/2 - 1))*...
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
        1/2*(-AlphaT_e/2 - 1)*((ti - abs(u(i)))/ti)^(-AlphaT_e/2 - 2)*(-1/ti*sign(u(i))*N(i, :)))*u(i)* ...
        (1 - GammaT_e*sin(atan(u_zeta(i)*(2/h)))^2)*(h/2) + ...
        ...
        N(i, :)'*2*MuT_e*(((ti - abs(u(i)))/ti)^(AlphaT_e - 1) + 1/2*((ti - abs(u(i)))/ti)^(-AlphaT_e/2 - 1))*N(i, :)*...
        (1 - GammaT_e*sin(atan(u_zeta(i)*(2/h)))^2)*(h/2) + ...
        ...
        N(i, :)'*2*MuT_e*(((ti - abs(u(i)))/ti)^(AlphaT_e - 1) + 1/2*((ti - abs(u(i)))/ti)^(-AlphaT_e/2 - 1))*u(i)*...
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

%% Defined auxiliary functions
function plot_result(ds, xs, PropertyTable)
% Plotting needle shape
% x = -sb:h:l;
figure('Name', sprintf('FEM Solution with Load Stepping'));
plt = plot(xs, rand(1, size(xs, 2)), 'k-');
grid on;
set(plt, 'YData', ds);
axis equal;

% Plotting patches
intervals = PropertyTable.Interval;
num_inter = size(intervals, 1);
Ps = [];
titles = [];
Py_Lims = get(gca, 'YLim');
for j = 1:num_inter
    Px = [intervals{j}(1), intervals{j}(2), intervals{j}(2), intervals{j}(1)];
    Py = [Py_Lims(1), Py_Lims(1), Py_Lims(2), Py_Lims(2)];
    if j == 1
        Ps = patch(Px, Py, 'w', 'FaceAlpha', 0.2);
    else
        Ps = [Ps, patch(Px, Py, j/num_inter, 'FaceAlpha', 0.2)];
    end
    titles = [titles; sprintf("Mu_{T%d} = %.0f Pa, Alpha_{T%d} = %.2f, Gamma_{T%d} = %.2f", ...
        j - 1, PropertyTable.MuT(j)*10^6, ...
        j - 1, PropertyTable.AlphaT(j), ...
        j - 1, PropertyTable.GammaT(j))];
end
legend(Ps, titles);
axis tight
end