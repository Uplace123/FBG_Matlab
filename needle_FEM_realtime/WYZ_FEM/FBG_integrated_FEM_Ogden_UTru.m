function [ds, ks, xs] = FBG_integrated_FEM_Ogden_UTru(sb, l, Db, Kb, ti, ...
    Nel, Mu, Alpha, Interval, curvatures, AA_lcn, d_init)
%% Initialization
% sb = 0;
% l = 165;
% Db = 0;
% Kb = 0;
% ti = 25;
% Nel = sb + l;
% Mu = 0;
% Alpha = 5;

% Add Cholesky inverse

% Constants
L = sb + l;
FEM_params;
% Kb = -Db/(sb + rcm);

% Interval = {[-sb, 0]; [0, l]};
MuT = Mu*10^-6; % Pa, but in mm^2; 1Pa = 1e-6 N/mm^2; 1kPa = 1e-3 N/mm^2
AlphaT = Alpha; % need abs(alpha) > 1
GammaT = zeros(size(MuT));
PropertyTable = table(Interval, MuT, AlphaT, GammaT);
% ti = 25; % Initial length of undeformed tissue. Will affect solution

% FEM-specific constants
% Nel = 10; % Total umber of elements
Nen = 4; % Number of element DOF
DOF = 1:(2*Nel + 2);
nDOF = length(DOF);
d = d_init; % Initial guess of solution
h = L/Nel; % Finite element size
EBC_idx = [1; 2]; % Displacement and slope of first element left node is prescribed
freeDOF = DOF;
freeDOF(EBC_idx) = [];

% Construction of LM
LM = zeros(Nen, Nel);
for e = 1:Nel
    LM(1, e) = 2*e - 1;
    LM(2, e) = 2*e - 0;
    LM(3, e) = 2*e + 1;
    LM(4, e) = 2*e + 2;
end

% Curvature inputs from FBG
AA_crv = 1e-3 * curvatures'; % curvature at AAs in mm
AA_er = round(AA_lcn./h) + 1; % elements where the left-moment is fixed
AA_crv = AA_crv(AA_er >= 0); % curvatures that have negative element indices are skipped
AA_er = AA_er(AA_er >= 0); % elements that have negative indices are skipped

%% FEM Main
K = zeros(nDOF, nDOF);
F = zeros(nDOF, 1);
P = zeros(nDOF, 1);
% Apply EBC
d(EBC_idx) = [Db; Kb];
% Loop over each element
for e = 1:Nel
    % Get local nodal values from global d
    d_i_local = d(2*e - 1:2*e + 2);
    [ke, pe] = compute_element_matrix(d_i_local, ti, E, I, ...
                PropertyTable, e, h, -sb, ...
                AA_er, AA_crv);
    % Global assembly process
    % Using vectors instead of FOR loops
    K(LM(1:4, e), LM(1:4, e)) = K(LM(1:4, e), LM(1:4, e)) + ke;
    P(LM(1:4, e)) = P(LM(1:4, e)) + pe;
end

% Newton's method
% Use only free DOF from the list of DOF to compute d
delta_d = invChol_mex(K(freeDOF, freeDOF))*(F(freeDOF, 1)-P(freeDOF, 1));
d(freeDOF, 1) = d(freeDOF, 1) + delta_d;

% output only the displacement
ds = extract_dist_from_d(d);
ks = extract_slop_from_d(d);
dds = diff(ds);
xs = zeros(1, Nel + 1);
xs(1) = -sb;
for i = 2:(Nel + 1)
    xs(i) = xs(i - 1) + sqrt(h^2 - dds(i - 1)^2);
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
    PropertyTable, e, h, x_begin, ...
    AA_er, AA_crv)
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
% modify element internal force vector based on FBG
FBG_idx_r = find(e == AA_er);
if FBG_idx_r
    pe(2) = AA_crv(FBG_idx_r)*E*I;
end

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
