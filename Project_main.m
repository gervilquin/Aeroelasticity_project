%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aeroelastic model of a Pazy wing
% Authors: Gerard Villalta & Antoni Alberto
% Date: 18/11/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc, clear, close all
addpath(genpath('functions'))

%% 1. Data input
% Geometrical data
c = 0.1; % chord of the wing (m)

y_sec = [0.0; 38.4; 76.8; 115.2; 153.6; 192.0; 230.4; 268.8; 307.2; 345.6; ...
         384.0; 422.4; 460.8; 499.2; 537.6; 541.6; 545.6; 550.0]*1e-3; % Section Y coordinate (m)
Nsec = length(y_sec)-1;             % Number of sections
Nesec = 5;                          % Number of elements per section
RibExist = [0, ones(1,Nsec-1), 0];  % Logical matrix for sections that contain ribs

% Material properties
material.Al.rho = 2795;       % Aluminium density (kg/m^3)
material.Al.E = 71000*1e6;    % Aluminium Young's modulus (Pa)
material.Al.v = 0.33;         % Aluminium Poisson coefficient (-)
material.Nylon.rho = 930;     % Nylon density (kg/m^3)
material.Nylon.E = 1700*1e6;  % Nylin Young's modulus (Pa)
material.Nylon.v = 0.394;     % Nylon Poisson coefficient (-)

% Structural properties (Pseudo-Experimental) %%%%% REVIEW VALUES %%%%%
x_sc = 0.43;           % Position of shear center (%chord)
EI = 5.04517;%6500*1e-3; % Mean value of flexural rigidity (N·m^2)
GJ = 6.5977;%5500*1e-3;  % Mean value of torsional rigidity (N·m^2/rad)
solve_static = true;   % Logical input to active/deactivate static solver

% Aerodynamic properties
x_ac = 1/4;  % Position of aerodynamic center (%chord)
x_col = 3/4; % Position of collocation point (%chord)
t = 18;      % Thickness of the airfoil (%chord)
U_inf = 60;  % Freestream velocity (m/s)
AoA = 10;     % Wing angle of attack (º)
rho_inf = 1; % Reference air density (kg/m^3)

%% 2. Mesh construction
% Mesh definition
y_el = ComputeYcoordinates(y_sec,Nsec,Nesec);  % Elements' nodes Y coordinate (m)
Nnodes = length(y_el);
Nel = Nnodes - 1;

% Connectivity matrices definition
Tn = ConnectivityElements(Nel);
Ts = ConnectivitySubsets(Tn,Nsec,Nesec);

% Dirichlet boundary conditions
Up = [  0 1 1;
        0 1 2;
        0 1 3];

%% 3. Structural modelling
% Define stiffness matrix
K = ComputeKmatrix(y_el,Tn,EI,GJ);
% Define mass matrix
M = ComputeMmatrix(y_el,Tn,Ts,material,t,c,x_sc,RibExist);

%% 4. Aerodynamics modelling
% Compute aerodynamic center's and colocation's points 
col_pos = ComputePointCoordinatesAerodynamics(y_el,c,x_col);
vortex_coord = ComputeVortexLinesNodes(y_el,x_ac,c);

% Compute Influence matrix
A = ComputeAmatrix(col_pos,vortex_coord);
% Compute element's size matrix
S = ComputeElementSize(y_el);

%% 5. Aeroelastic linear coupling
% Compute the aerodynamic coupling matrices
[I_au_1, I_au_2, I_au_3] = ComputeDisplacementsCoupling(y_el,Tn,U_inf);
% Compute the structural coupling matrix
e = (x_ac - x_sc)*c;
I_fL = ComputeForcesCoupling(y_el,Tn,e);

%% 6. Static case solver
% Compute the indeces of the free fix DOF
[Ip,If,u_static] = ComputeBoundaryConditions(Up,3*Nnodes);
% Compute force vector
[F,L] = ComputeFvector(Nel,U_inf,rho_inf,AoA,S,A,I_fL);

if solve_static == true
    [u] = StaticSolver(K,F,u_static,If,Ip);
    plotStaticSolution(u,F,L,y_el,U_inf,rho_inf,AoA,true)
end

%% 7. Divergence solver
Uinf_ = linspace(0.1,150,100);
U_diverg = [0];
w_tip_Uinf = [0];


for i=1:length(Uinf_)
    U_inf = Uinf_(i);

    [Meff,Ceff,Keff] = ComputeEffectiveMatrix(M,K,U_inf,Parameters_list);

 
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
Nm = 6; %first 6 modes

% Obtain the first eigenvectos and eigenvalues
[V,D] = eigs(Kf,Mf,Nm,'sr');

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
    plot(y_nodal,Phi(2:3:end,i),'DisplayName',string(round(w2(i),0)))
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

%Uinf_ = logspace(-10,-1,100);
Uinf_ = linspace(0.1,U_diverg(end),3);
%Uinf_ = [1];
p_values = zeros(length(Uinf_),1);
p_values_red = zeros(length(Uinf_),1);
N_modes = 30;
p_values_collect = zeros(length(Uinf_),N_modes);

% Get the eigenvalues of M and K
[eig_vector, eig_value] = eigs(K,M,N_modes,'sm');
eig_vector = eig_vector(If,:);


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

    [eig_vector_p_red, eig_value_p_red] = eigs(D,N_modes,'sm');



    p_values_red(i) = max(real(-1./diag(eig_value_p_red)));
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





