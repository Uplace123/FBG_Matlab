% generate FEM constant parameters as mat file

% constant used in FEM
sb = 0; % l1
l = 165; % l2
Db = 0;
Kb = 0; % Kb = -Db/(sb + rcm);
ti = 25; % ti = 25; % Initial length of undeformed tissue. Will affect solution
Mu = 0;
Alpha = 5;
E = 250*1000; % 200GPa but in mm^2
OD = 1.27;
I = pi/4*(OD/2)^4; % in mm^4
Interval = {[-sb, 0]; [0, l]};
L = sb + l; % Constants whole length
MuT = [0; Mu]*10^-6; % Pa, but in mm^2; 1Pa = 1e-6 N/mm^2; 1kPa = 1e-3 N/mm^2
AlphaT = [0; Alpha]; % need abs(alpha) > 1 
GammaT = [0; 0];
inner_iter = 0; % Initialize inner_iter counter
outer_iter = 0; % Initialize outer_iter counter
converged = 0;
max_inner_iter = 5; % Max number of iterations for Newton's method
max_outer_iter = 50; % Max number of iterations for load stepping
tol = 1e-3; % Convergence criterion
load_ratio = 1;
Nen = 4; % Number of element DOF
PropertyTable = table(Interval, MuT, AlphaT, GammaT);
EBC_cur = zeros(2, 1); % For load stepping
EBC_converged = zeros(2, 1); % For load stepping from previously converged EBC
Nel = sb + l; % total numbr of element
DOF = 1:(2*Nel + 2);
freeDOF = DOF;
h = L/Nel; % Finite element size
nDOF = length(DOF);
d = zeros(nDOF, 1); % Initial guess of solution
EBC = [1; 2]; % Displacement and slope of first element left node is prescribed
freeDOF(EBC) = [];
LM = zeros(Nen, Nel); % Construction of LM
for e = 1:Nel
    LM(1, e) = 2*e - 1;
    LM(2, e) = 2*e - 0;
    LM(3, e) = 2*e + 1;
    LM(4, e) = 2*e + 2;
end

%AA_lcn = [65, 100, 135,155]; % location of AA wrt to base
AA_lcn = [65,100,135];
save("FEM_params.mat");