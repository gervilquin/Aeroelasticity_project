%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aeroelastic model of a Pazy wing
% Authors: Gerard Villalta & Antoni Alberto
% Date: 18/11/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc, clear, close all
addpath('functions\')
%% User inputs

    % Aerodynamic
    U_inf = 60; %m/s
    AoA = 7;
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
y_sections(:,2) = y_sections(:,2)*1e-3; % Coordiantes in meters

Nsec = length(y_sections(:,2))-1; % number of panels

Nel = Nsec*Neset; % total number of elements

% Centers >Previously computed<
x_sc = 0.43;
x_ac = 1/4;
x_col = 3/4;

% Geometry
chord = 100*1e-3; %m

% Material properties >Previously computed<
EI = 5.04517;%6500*1e-3;
GJ = 6.5977;%5500*1e-3;

%% 2. Structural modelling

% Material properties
material.Al.rho = 2795;
material.Al.E = 71000*1e6;
material.Al.v = 0.33;
material.Nylon.rho = 930;
material.Nylon.E = 1700*1e6;
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
M = def_M_matrix(y_nodal,Neset,Nsec,Nnod,Tn,Ts,material,naca_2dt,chord,x_sc*chord);


% plot matrices
figure()
title("K")
x_surf = linspace(1,length(K),length(K));
surf(x_surf,x_surf,K,'EdgeColor','none')
%set(gca,'zscale','log')

figure()
title("M")
x_surf = linspace(1,length(M),length(M));
surf(x_surf,x_surf,M,'EdgeColor','none')
%set(gca,'zscale','log')

%% 3. Aerodynamics modelling

% Define the chord vector
chord_y = ones(1,length(y_nodal))*chord;

% Compute area of the elements
%S = compute_element_surface(chord_y,y_nodal);
S = compute_element_length(y_nodal);

% Compute chord/4 point and the colocation point spanwise location
ac_pos = compute_aero_point(y_nodal,chord_y,x_ac);%*1e-3;
col_pos = compute_aero_point(y_nodal,chord_y,x_col);%*1e-3;
segment_coor = compute_segment_coordinate(y_nodal,ac_pos,chord_y);%*1e-3;

% Compute Influence matrix
A_aero = compute_A_matrix(col_pos,segment_coor);

%% 4. Aeroelastic linear coupling

% Compute the aerodynamic coupling matrices
[I_au_0,I_au_1,I_au_2] = compute_I_au(U_inf,Nnod,Nel,Tn);

% Compute the structural coupling matrix
I_fL = compute_I_fL(Nnod,Nel,-(x_ac-x_sc)*chord);

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

I_au_0_f = I_au_0(:,If);
I_au_1_f = I_au_1(:,If);
I_au_2_f = I_au_2(:,If);

I_fL_f = I_fL(If,:);

% Collect all parameters in a list
Parameters_list.S = S;
Parameters_list.rho = rho_inf;
Parameters_list.Nnod = Nnod;
Parameters_list.Nel = Nel;
Parameters_list.Tn = Tn;
Parameters_list.If = If;
Parameters_list.I_fl = I_fL;
Parameters_list.A_aero = A_aero;

%% 5 Static case solution
if solve_static == true
    [Meff,Ceff,Keff] = compute_efective_matrices(M,K,U_inf,Parameters_list);
    solve_static_case(Nnod,y_nodal,u_static,If,Ip,I_fL_f,S,A_aero,Keff,U_inf,rho_inf,AoA,true);
    %solve_static_case(Nnod,y_nodal,u_static,If,Ip,I_fL_f,S,A_aero,K,U_inf,rho_inf,AoA,true);
end

%% Compute divergence
Uinf_ = linspace(0.1,150,100);
U_diverg = [0];
w_tip_Uinf = [0];


for i=1:length(Uinf_)
    U_inf = Uinf_(i);

    [Meff,Ceff,Keff] = compute_efective_matrices(M,K,U_inf,Parameters_list);

 
    %close all
    [theta,w,gamma] = solve_static_case(Nnod,y_nodal,u_static,If,Ip,I_fL_f,S,A_aero,Keff,U_inf,rho_inf,AoA,false);
    
    if (w(end-1) - w_tip_Uinf(end) > 0)

        w_tip_Uinf(end+1) = w(end-1);
        U_diverg(end+1) = U_inf;
    end

end

figure()
hold on
plot(U_diverg,w_tip_Uinf)
xlabel("$U_{\infty}$","Interpreter","latex")
ylabel("$w^{tip}$","Interpreter","latex")
grid on
grid minor
legend([strcat("Divergence $U_{\infty}$ = ",string(U_diverg(end))," m/s")],'location','northwest',"Interpreter","latex")
hold off

%% Modal analysis

% define the number of modes that want to be returned
Nm = 4; %first 6 modes

% Obtain the first eigenvectos and eigenvalues
[V,D] = eigs(Kf,Mf,Nm,'sm');

% Obtain the natural frequencies and the vibration modes
Phi = zeros(Nnod*3,Nm); 
w2 = zeros(1,Nm);

for k =1:length(V(1,:))
    Phi(If,k) = V(:,k)/sqrt(V(:,k)'*M(If,If)*V(:,k));
    w2(k) = D(k,k);
end

% plot modes
figure()

subplot(3,1,1)
hold on
for i =1:Nm
    plot(y_nodal,Phi(1:3:end,i),'DisplayName',string(round(w2(i),0)))
end
%legend()
hold off
ylabel("Twist","Interpreter","latex")
subplot(3,1,2)
hold on
for i =1:Nm
    plot(y_nodal,Phi(2:3:end,i),'DisplayName',string(round(sqrt(w2(i))/2/pi,0)))
end
legend(Location='bestoutside')
hold off
ylabel("U vertical","Interpreter","latex")
subplot(3,1,3)
hold on
for i =1:Nm
    plot(y_nodal,Phi(3:3:end,i),'DisplayName',string(round(w2(i),0)))
end
%legend()
hold off
ylabel("Deflection","Interpreter","latex")
xlabel("Y [mm]",'Interpreter','latex')




%% 5. Aeroelastic solver

Uinf_ = linspace(0.1,0.2*U_diverg(end),200);

%Uinf_ = [1];
p_values = zeros(length(Uinf_),1);
p_values_red = zeros(length(Uinf_),1);
N_modes = 6;
p_values_collect = zeros(length(Uinf_),N_modes);

% Get the eigenvalues of M and K
[eig_vector, eig_value] = eigs(Kf,Mf,N_modes,'sm');
%eig_vector = eig_vector(If,:);


for i = 1:length(Uinf_)
    U_inf = Uinf_(i);
    rho_inf = 1;

%     % Compute Aerodynamic matrices
%     S_aero = -U_inf^2*rho_inf*S;
% 
%     % Compute coupling matrices
%     [I_au_0,I_au_1,I_au_2] = compute_I_au(U_inf,Nnod,Nel,Tn);
%     %I_fL = compute_I_fL(Nnod,Nel,x_sc-x_ac);
% 
%     I_au_0 = I_au_0(:,If);
%     I_au_1 = I_au_1(:,If);
%     I_au_2 = I_au_2(:,If);
% 
%     % Compute aero mass, stiffness and damping matrices
%     M_a = I_fL*(S_aero*inv(A_aero))*I_au_2;
%     C_a = I_fL*(S_aero*inv(A_aero))*I_au_1;
%     K_a = I_fL*(S_aero*inv(A_aero))*I_au_0;
%     
%     % Compute efective matrices
%     Meff = Mf + M_a; 
%     Ceff = C_a;
%     Keff = Kf + K_a;

    [Meff,Ceff,Keff] = compute_efective_matrices(M,K,U_inf,Parameters_list);

    Meff = Meff(If,If);
    Ceff = Ceff(If,If);
    Keff = Keff(If,If);

%     % >> P method as a quadratic eigenvalue problem <<
%     B = [Ceff Meff; -1*eye(size(Keff)) zeros(size(Keff))];
%     A = [Keff zeros(size(Keff)); zeros(size(Keff)) eye(size(Keff))];
% 
%     [eig_vector_p, eig_value_p] = eigs(A,B,30,'sm');
% 
%     p_values(i) = max(real(-1./diag(eig_value_p)));

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

    [eig_vector_p_red, eig_value_p_red] = eigs(D,N_modes,'lm');
    %p_values_red(i) = max(real(1./diag(eig_value_p_red)));

    p_values_red(i) = max(real(1./diag(eig_value_p_red)));
    p_values_collect(i,:) = -1./diag(eig_value_p_red);
end


%% Plots
% figure()
% plot(Uinf_,p_values)
% grid on
% grid minor
% ylabel("$max(Re(p_i))$",'Interpreter','latex')
% xlabel("$U_{\infty}$",'Interpreter','latex')
% hold off

figure()
hold on
for i = 1:length(Uinf_)
    plot(real(p_values_collect(i,:)),imag(p_values_collect(i,:)))
end
hold off

figure()
plot(Uinf_,p_values_red)
grid on
grid minor
ylabel("$max(Re(p_i))$",'Interpreter','latex')
xlabel("$U_{\infty}$",'Interpreter','latex')
hold off

    





