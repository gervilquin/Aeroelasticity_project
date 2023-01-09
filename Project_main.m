%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aeroelastic model of a Pazy wing
% Authors: Gerard Villalta & Antoni Alberto
% Date: 18/11/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc, clear, close all
addpath('functions\')
%% User inputs

    % Aerodynamic
    U_inf = 1; %m/s
    AoA = 1;
    rho_inf = 1;

    % Number of elements between ribs
    Neset = 5;

    % To solve the static case
    solve_static = true;



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

Nel = Nsec*Neset; % total number of elements

% Centers >Previously computed<
x_sc = 0.43;
x_ac = 1/4;
x_col = 3/4;

% Geometry
chord = 100; %mm

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
Ts = ConnectivitySubsets(Tn,Nsec-1,Neset);

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
ac_pos = compute_aero_point(y_nodal,chord_y,x_ac)*1e-3;
col_pos = compute_aero_point(y_nodal,chord_y,x_col)*1e-3;
segment_coor = compute_segment_coordinate(y_nodal,ac_pos,chord_y)*1e-3;

% Compute Influence matrix
A_aero = compute_A_matrix(col_pos,segment_coor);

%% 4. Aeroelastic linear coupling

% Compute the aerodynamic coupling matrices
[I_au_0,I_au_1,I_au_2] = compute_I_au(U_inf,Nnod,Nel,Tn);

% Compute the structural coupling matrix
I_fL = compute_I_fL(Nnod,Nel,x_ac-x_sc);

%% Boundary conditions

% Fixed DOF
Up = [  0 1 1;
        0 1 2;
        0 1 3];

% Compute the indeces of the free fix DOF
[Ip,If,u_static] = compute_boundary_conditions(Up,3*Nnod);

% Split the matrices in free and fix
Kf = K(If,If);
Mf = M(If,If);

I_au_0 = I_au_0(:,If);
I_au_1 = I_au_1(:,If);
I_au_2 = I_au_2(:,If);

I_fL = I_fL(If,:);

%% 5 Static case solution
if solve_static == true
    solve_static_case(Nnod,y_nodal,u_static,If,Ip,I_fL,S,A_aero,K,U_inf,rho_inf,AoA,true);
end

%% Compute divergence
Uinf_ = linspace(0.1,300,100);
U_diverg = [];
w_tip_Uinf = [];

for i=1:length(Uinf_)
    U_inf = Uinf_(i);

    [I_au_0,I_au_1,I_au_2] = compute_I_au(U_inf,Nnod,Nel,Tn);
    I_au_0 = I_au_0(:,If);
    I_au_1 = I_au_1(:,If);
    I_au_2 = I_au_2(:,If);

    try 
        close all
        [theta,w,gamma] = solve_static_case(Nnod,y_nodal,u_static,If,Ip,I_fL,S,A_aero,K,U_inf,rho_inf,AoA,true);
        w_tip_Uinf(end+1) = w(end-1);
        U_diverg(end+1) = U_inf;
    catch
        pass
    end
end

figure()
plot(U_diverg,w_tip_Uinf)


%% 5. Aeroelastic solver

%Uinf_ = logspace(-10,-1,100);
Uinf_ = linspace(0.1,100,100);
p_values = zeros(length(Uinf_),1);
p_values_red = zeros(length(Uinf_),1);

% Get the eigenvalues of M and K
[eig_vector, eig_value] = eigs(Kf,Mf,10,'sm');

for i = 1:length(Uinf_)
    U_inf = Uinf_(i);
    rho_inf = 1;

    % Compute Aerodynamic matrices
    S_aero = -U_inf^2*rho_inf*S;

    % Compute coupling matrices
    [I_au_0,I_au_1,I_au_2] = compute_I_au(U_inf,Nnod,Nel,Tn);
    %I_fL = compute_I_fL(Nnod,Nel,x_sc-x_ac);

    I_au_0 = I_au_0(:,If);
    I_au_1 = I_au_1(:,If);
    I_au_2 = I_au_2(:,If);

    % Compute aero mass, stiffness and damping matrices
    M_a = I_fL*(S_aero*inv(A_aero))*I_au_2;
    C_a = I_fL*(S_aero*inv(A_aero))*I_au_1;
    K_a = I_fL*(S_aero*inv(A_aero))*I_au_0;
    
    % Compute efective matrices
    Meff = Mf + M_a; 
    Ceff = C_a;
    Keff = Kf + K_a;

    % >> P method as a quadratic eigenvalue problem <<
    B = [Ceff Meff; -1*eye(size(Keff)) zeros(size(Keff))];
    A = [Keff zeros(size(Keff)); zeros(size(Keff)) eye(size(Keff))];

    [eig_vector_p, eig_value_p] = eigs(A,B,30,'sm');

    p_values(i) = max(real(-1./diag(eig_value_p)));

    % >> reduced set
    Meff_red = eig_vector'*Meff*eig_vector;
    Ceff_red = eig_vector'*Ceff*eig_vector;
    Keff_red = eig_vector'*Keff*eig_vector;

%     B_red = [Ceff_red  Meff_red ; -1*eye(size(Keff_red )) zeros(size(Keff_red ))];
%     A_red = [Keff_red  zeros(size(Keff_red )); zeros(size(Keff_red )) eye(size(Keff_red))];
%     
%     [eig_vector_p_red, eig_value_p_red] = eigs(A_red,B_red,30,'sm');
    
    D = [Keff_red\Ceff_red Keff_red\Meff_red;
        -1*eye(size(Keff_red)) zeros(size(Keff_red))];

    [eig_vector_p_red, eig_value_p_red] = eigs(D,30,'sm');

    p_values_red(i) = max(real(-1./diag(eig_value_p_red)));
end


%% Plots
figure()
plot(Uinf_,p_values)
grid on
grid minor
ylabel("$max(Re(p_i))$",'Interpreter','latex')
xlabel("$U_{\infty}$",'Interpreter','latex')
hold off

figure()
plot(Uinf_,p_values_red)
grid on
grid minor
ylabel("$max(Re(p_i))$",'Interpreter','latex')
xlabel("$U_{\infty}$",'Interpreter','latex')
hold off





