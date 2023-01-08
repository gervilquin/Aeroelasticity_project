%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aeroelastic model of a Pazy wing
% Authors: Gerard Villalta & Antoni Alberto
% Date: 18/11/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc, clear, close all
addpath('functions\')

%% 1. Initialisation 

% Definition of the location of the ribs
%                   Rib     Node j      yj  
y_sections =       [0,       1,           0.0
                    1,       2,           38.4
                    2,       3,           76.8
                    3,       4,           115.2
                    4,       5,           153.6
                    5,       6,           192.0
                    6,       7,           230.4
                    7,       8,           268.8
                    8,       9,           307.2
                    9,       10,          345.6
                    10,      11,          384.0
                    11,      12,          422.4
                    12,      13,          460.8
                    13,      14,          499.2
                    14,      15,          537.6
                    15,      16,          541.6
                    16,      17,          545.6
                    0,       18,          550.0 ];

y_sections = y_sections(:,2:3);

Nsec = length(y_sections(:,2))-1; % number of panels
Neset = 5;%<<<<<<<<<< Input 

Nel = Nsec*Neset; % total number of elements

% Centers >Previously computed<
x_sc = 0.43;
x_ac = 1/4;
x_col = 3/4;

% Geometry
chord = 100; %mm

% Aerodynamic
U_inf = 100; %m/s

% Material properties >Previously computed<
EI = 6500;
GJ = 5500;

%% 2. Structural modelling

% Material properties
material.Al.rho = 2795;
material.Al.E = 71000;
material.Al.v = 0.33;
material.Nylon.rho = 930;
material.Nylon.E = 1700;
material.Nylon.v = 0.394;
naca_2dt = 18;

% Define the nodal coordinates
y_nodal = compute_x_nodal(y_sections(:,2),Neset);
Nnod = length(y_nodal); % number

% define the connectivity matrixs
Tn = ConnectivityElements(Nnod-1);
Ts = ConnectivitySubsets(Tn,Nsec-1,Neset); % PREGUNTA -> Perquè Ts(end) != Tn(end)

% Define stiffness matrix
K = def_K_matrix(y_nodal,Neset,Nsec,Nnod,Tn,EI,GJ);


% Define mass matrix
M = def_M_matrix(y_nodal,Neset,Nsec,Nnod,Tn,Ts,material,naca_2dt,chord,x_sc);


%% 3. Aerodynamics modelling

% Define the chord vector
chord_y = ones(1,length(y_nodal))*chord;

% Compute area of the elements
%S = compute_element_surface(chord_y,y_nodal);
S = compute_element_length(y_nodal);

% Compute chord/4 point and the colocation point spanwise location
ac_pos = compute_aero_point(y_nodal,chord_y,x_ac);
col_pos = compute_aero_point(y_nodal,chord_y,x_col);
segment_coor = compute_segment_coordinate(y_nodal,ac_pos,chord_y);

% Compute Influence matrix
A_aero = compute_A_matrix(col_pos,segment_coor);

%% 4. Aeroelastic linear coupling

% Compute the aerodynamic coupling matrices
[I_au_0,I_au_1,I_au_2] = compute_I_au(U_inf,Nnod,Nel,Tn);

% Compute the structural coupling matrix
I_fL = compute_I_fL(Nnod,Nel,x_sc-x_ac);

% u = zeros(3*Nnod,1);
% u_dot = zeros(3*Nnod,1);
% u_dotdot = zeros(3*Nnod,1);
% 
% alpha = I_au_0*u + I_au_1*u_dot + I_au_2*u_dotdot;
% L = -U_inf^2*S\A_aero*alpha;
% F_nod = I_fL*L;

%% Static case solution
solve_static = false;
if solve_static == true
    solve_static_case(Nnod,y_nodal,I_au_0,I_au_1,I_au_2,I_fL,S,A_aero,K,U_inf)
end


%% 5. Aeroelastic solver

Uinf_ = logspace(-10,1,100);
p_values = zeros(length(Uinf_),1);

% Get the eigenvalues of M and K
%[eig_vector, eig_value] = eigs(K,M,20,'sm');

for i = 1:length(Uinf_)
    U_inf = Uinf_(i);
    rho_inf = 1;

    % Compute Aerodynamic matrices
    S_aero = -U_inf^2*rho_inf*S;

    % Compute coupling matrices
    [I_au_0,I_au_1,I_au_2] = compute_I_au(U_inf,Nnod,Nel,Tn);
    I_fL = compute_I_fL(Nnod,Nel,x_sc-x_ac);

    % Compute aero mass, stiffness and damping matrices
    M_a = I_fL*(S_aero*inv(A_aero))*I_au_2;
    C_a = I_fL*(S_aero*inv(A_aero))*I_au_1;
    K_a = I_fL*(S_aero*inv(A_aero))*I_au_0;
    
    % Compute efective matrices
    Meff = M + M_a; 
    Ceff = C_a;
    Keff = K + K_a;

    % >> P method as a quadratic eigenvalue problem <<
    B = [Ceff Meff; -1*eye(size(Keff)) zeros(size(Keff))];
    A = [K zeros(size(Keff)); zeros(size(Keff)) eye(size(Keff))];

    [eig_vector, eig_value] = eigs(A,B,30,'sm');

    p_values(i) = max(real(-1./diag(eig_value)));
    
end


%% Plots
figure()
plot(Uinf_,p_values)
grid on
grid minor
ylabel("$max(Re(p_i))$",'Interpreter','latex')
xlabel("$U_{\infty}$",'Interpreter','latex')





