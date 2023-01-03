%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aeroelastic model of a Pazy wing
% Authors: Gerard Villalta & Antoni Alberto
% Date: 18/11/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc, clear, close all
addpath('functions\')

%% 1. Initialisation 

% Elmements definition
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

nsec = length(y_sections(:,2))-1; % number of panels
neset = 5;%<<<<<<<<<< Input 

nel = nsec*neset; % total number of elements

% Centers
x_sc = 0.43;
x_ac = 1/4;
x_col = 3/4;

% Geometry
chord = 100; %mm

% Aerodynamic
U_inf = 1; %m/s

% Material properties
EI = 6500;
GJ = 5500;

%% 2. Structural modelling

material.Al.rho = 2795;
material.Al.E = 71000;
material.Al.v = 0.33;
material.Nylon.rho = 930;
material.Nylon.E = 1700;
material.Nylon.v = 0.394;
naca_2dt = 18;

% Define the nodal coordinates
y_nodal = compute_x_nodal(y_sections(:,2),neset);
n = length(y_nodal);

% define the connectivity matrixs
Tn = ConnectivityElements(n-1);
Ts = ConnectivitySubsets(Tn,nsec-1,neset);

% Define stiffness matrix
K = sparse(3*n,3*n);
for i = 1:neset*nsec
    I = zeros(2*3,1); % (nnodes_el*nDOFs, 1)
    l = y_nodal(Tn(i,2)) - y_nodal(Tn(i,1));
    for k = 1:3 % 3 DOFs
        I(k,1) = 3*(Tn(i,1)-1)+k;
        I(k+3,1) = 3*(Tn(i,2)-1)+k;
    end
    Kt = (GJ/l)*[ 1 0 0 -1 0 0;
                  0 0 0  0 0 0;
                  0 0 0  0 0 0;
                 -1 0 0  1 0 0;
                  0 0 0  0 0 0;
                  0 0 0  0 0 0];
    Kb = (EI/l^3)*[0  0    0   0   0   0   ;
                   0  12  6*l  0  -12  6*l ; 
                   0 6*l 4*l^2 0 -6*l 2*l^2;
                   0  0    0   0   0   0   ;
                   0 -12 -6*l  0   12  -6*l;
                   0 6*l 2*l^2 0 -6*l 4*l^2];
    K(I,I) = K(I,I) + Kt + Kb;
end

% Define mass matrix
M = sparse(3*n,3*n);
for i = 1:neset*nsec
    I = zeros(2*3,1); % (nnodes_el*nDOFs, 1)
    l = y_nodal(Tn(i,2)) - y_nodal(Tn(i,1));
    for k = 1:3 % 3 DOFs
        I(k,1) = 3*(Tn(i,1)-1)+k;
        I(k+3,1) = 3*(Tn(i,2)-1)+k;
    end
    [rhoA] = compute_avg_density(material,naca_2dt,chord);
    [x_cm] = compute_mass_center(material,rhoA,naca_2dt,chord);
    [I_cm] = compute_inertia_mass_center(x_cm);

    Mt = (Icm*l/6)*[2 0 0 1 0 0;
                    0 0 0 0 0 0;
                    0 0 0 0 0 0;
                    1 0 0 2 0 0;
                    0 0 0 0 0 0;
                    0 0 0 0 0 0];
    Mb = (rhoA*l/420)*[0   0       0   0    0       0 ;
                       0  156    22*l  0   54    -13*l;
                       0  22*l  4*l^2  0  13*l  -3*l^2;
                       0   0       0   0    0       0 ;
                       0   54    13*l  0   156   -22*l;
                       0 -13*l -3*l^2  0 -22*l   4*l^2];
    d = x_cm - x_sc;
    d_m = [ 1  0  0  0  0  0;
         -d  1  0  0  0  0;
          0  0  1  0  0  0;
          0  0  0  1  0  0;
          0  0  0 -d  1  0;
          0  0  0  0  0  1];
    M(I,I) = M(I,I) + d_m'*(Mt + Mb)*d;
end


%% 3. Aerodynamics modelling



% Define the chord vector
chord_y = ones(1,length(y_nodal))*chord;

% Compute area of the elements
S_aero = compute_element_surface(chord_y,y_nodal);

% Compute chord/4 point and the colocation point span wise location
ac_pos = compute_aero_point(y_nodal,chord_y,x_ac);
col_pos = compute_aero_point(y_nodal,chord_y,x_col);
segment_coor = compute_segment_coordinate(y_nodal,ac_pos,chord_y);

% Compute Influence matrix
A_aero = compute_A_matrix(col_pos,segment_coor);

% % Compute Lift
% alpha_test = 1*pi/180;
% alpha_array = alpha_test*ones(length(S_),1);
% 
% Lift = -1*U_inf^2*diag(S_/10^6)*inv(A_aero)*alpha_array;
% 
% figure()
% plot(ac_pos(:,2),Lift)




%% 4. Aeroelastic linear coupling

[I_au_0,I_au_1,I_au_2] = compute_I_au(U_inf,n,nel,Tn);
%I_fL = compute_I_fL(x_sc-x_ac);

u = zeros(3*n,1);
u_dot = zeros(3*n,1);
u_dotdot = zeros(3*n,1);

alpha = I_au_0*u + I_au_1*u_dot + I_au_2*u_dotdot;
L = -U_inf^2*S_aero\A_aero*alpha;

%% 5. Aeroelastic solver



