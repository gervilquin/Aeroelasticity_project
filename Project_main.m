%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aeroelastic model of a Pazy wing
% Authors: Gerard Villalta & Antoni Alberto
% Date: 18/11/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc, clear, close all
addpath('functions\')

% Case solution active/deactivate
solve_static = false; 
solve_diverge = false;
solve_modal = false;
solve_flutter = true;

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


% Aerodynamic properties
x_ac = 1/4;  % Position of aerodynamic center (%chord)
x_col = 3/4; % Position of collocation point (%chord)
t = 18;      % Thickness of the airfoil (%chord)
U_inf = 60;  % Freestream velocity (m/s)
AoA = 10;     % Wing angle of attack (º)
rho_inf = 1.3; % Reference air density (kg/m^3)

%% 2. Structural modelling
% Structural mesh definition
y_el = ComputeYcoordinates(y_sec,Nsec,Nesec);  % Elements' nodes Y coordinate (m)
Nnodes = length(y_el);
Nel = Nnodes - 1;

% Connectivity matrices definition
Tn = ConnectivityElements(Nel);
Ts = ConnectivitySubsets(Tn,Nsec,Nesec);

% Define stiffness matrix
K = ComputeKmatrix(y_el,Tn,EI,GJ);
% Define mass matrix
M = ComputeMmatrix(y_el,Tn,Ts,material,t,c,x_sc*c,RibExist);

% % plot matrices
% figure()
% title("K")
% x_surf = linspace(1,length(K),length(K));
% surf(x_surf,x_surf,K,'EdgeColor','none')
% set(gca,'zscale','log')
% 
% figure()
% title("M")
% x_surf = linspace(1,length(M),length(M));
% surf(x_surf,x_surf,M,'EdgeColor','none')
% set(gca,'zscale','log')

%% 3. Aerodynamics modelling
% Compute aerodynamic center's and colocation's points 
ac_pos =  x_ac*c;
col_pos = ComputePointCoordinatesAerodynamics(y_el,c,x_col);
vortex_coord = ComputeVortexLinesNodes(y_el,ac_pos,c);

% Compute Influence matrix
A = ComputeAmatrix(col_pos,vortex_coord);
% Compute element's size matrix
S = ComputeElementSize(y_el);

%% 4. Aeroelastic linear coupling
% Compute the aerodynamic coupling matrices
[I_au_0, I_au_1, I_au_2] = ComputeDisplacementsCoupling(y_el,Tn,U_inf);
% Compute the structural coupling matrix
e = (x_ac - x_sc)*c;
I_fL = ComputeForcesCoupling(y_el,Tn,e);

%% Boundary conditions
% Fixed DOF
Up = [  0 1 1;
        0 1 2;
        0 1 3];

% Compute the indeces of the free fix DOF
[Ip,If,u_static] = ComputeBoundaryConditions(Up,3*Nnodes);

% % Split the matrices in free and fix %%%% TO DO JUST BEFORE EIGS? %%%%
% Kf = K(If,If);
% Mf = M(If,If);
% I_au_1_f = I_au_1(:,If);
% I_au_2_f = I_au_2(:,If);
% I_au_3_f = I_au_3(:,If);
% I_fL_f = I_fL(If,:);

%% 5 Static case solution
if solve_static == true
    [Meff,Ceff,Keff] = ComputeEffectiveMatrix(K,M,y_el,Tn,U_inf,rho_inf,A,S,I_fL);
    StaticSolver(y_el,u_static,If,Ip,I_fL,S,A,Keff,U_inf,rho_inf,AoA,true);
end

%% 5.1 Static solution with aerocoupling
if solve_static == true
    [Meff,Ceff,Keff] = ComputeEffectiveMatrix(K,M,y_el,Tn,U_inf,rho_inf,A,S,I_fL);
    AeroStaticSolver(y_el,u_static,If,Ip,I_fL,I_au_0,S,A,Keff,U_inf,rho_inf,AoA,true);
end

%% Compute divergence
if solve_diverge == true
    Uinf_ = linspace(0.1,150,100);
    U_diverg = [0];
    w_tip_Uinf = [0];
    
    
    for i=1:length(Uinf_)
        U_inf = Uinf_(i);
    
        [Meff,Ceff,Keff] = ComputeEffectiveMatrix(K,M,y_el,Tn,U_inf,rho_inf,A,S,I_fL);
    
        [theta,w,gamma] = AeroStaticSolver(y_el,u_static,If,Ip,I_fL,I_au_0,S,A,Keff,U_inf,rho_inf,AoA,false);
        
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
end

%% Modal analysis
if solve_modal == true
    % define the number of modes that want to be returned
    Nm = 6; %first 6 modes
    
    % Obtain the first eigenvectos and eigenvalues
    [V,D] = eigs(K(If,If),M(If,If),Nm,'sm');
    
    % Obtain the natural frequencies and the vibration modes
    Phi = zeros(Nnodes*3,Nm); 
    w2 = zeros(1,Nm);
    freq = zeros(1,Nm); 
    
    for k =1:length(V(1,:))
        Phi(If,k) = V(:,k)/sqrt(V(:,k)'*M(If,If)*V(:,k));
        w2(k) = D(k,k);
        freq(k) = sqrt(w2(k))/2/pi;
    end
    
    % plot modes
    figure()
    tcl = tiledlayout(3,1);
    
    % Plot mode shpaes in the twist
    nexttile(tcl)
    hold on
    for i =1:Nm
        plot(y_el,Phi(1:3:end,i))
    end
    
    grid on
    grid minor
    hold off
    ylabel("Twist","Interpreter","latex")
    
    
    %Plot mode shapes in the vertical displacement 
    nexttile(tcl)
    hold on
    for i =1:Nm
        plot(y_el,Phi(2:3:end,i))
    end
    grid on
    grid minor
    hold off
    ylabel("U vertical","Interpreter","latex")

    %Plot the mode shapes in deflection
    nexttile(tcl)
    hold on
    for i =1:Nm
        plot(y_el,Phi(3:3:end,i),'DisplayName',strcat("Mode ",string(i),": ",string(round(freq(i),0))," Hz"))
    end
    grid on
    grid minor
    hold off
    ylabel("Deflection","Interpreter","latex")
    xlabel("Y [mm]",'Interpreter','latex')

    % Add legend to the side
    Lgnd = legend('show');
    Lgnd.Layout.Tile = 'East';

end



%% 5. Flutter solver
if solve_flutter == true
%Uinf_ = logspace(-10,-1,100);
%Uinf_ = linspace(0.1,U_diverg(end),3);
Uinf_ = linspace(0.1,100,100);

% Get the eigenvalues of M and K
N_reduced = 30;
N_modes = 10;
[Vr, Dr] = eigs(K(If,If),M(If,If),N_reduced,'sm');

% Initialize matrices
p_values = zeros(length(Uinf_),1);
p_values_collect = zeros(length(Uinf_),N_modes);

for i = 1:length(Uinf_)
    U_inf = Uinf_(i);

    % Compute efective matrix
    [Meff,Ceff,Keff] = ComputeEffectiveMatrix(K,M,y_el,Tn,U_inf,rho_inf,A,S,I_fL);

    % Select the Free DOFs
    Meff = Meff(If,If);
    Ceff = Ceff(If,If);
    Keff = Keff(If,If);

    % Order reduction
    Meff_red = Vr'*Meff*Vr;
    Ceff_red = Vr'*Ceff*Vr;
    Keff_red = Vr'*Keff*Vr;

    % Compute D matrix
    D = [Keff_red\Ceff_red Keff_red\Meff_red;
        -1*eye(size(Keff_red)) zeros(size(Keff_red))];

    % Compute eigen values
    [Vd, Dd] = eigs(D,N_modes,'sm');

    p_values(i) = max(real(-1./diag(Dd)));
    p_values_collect(i,:) = -1./diag(Dd);
end


%% Plots

% figure()
% hold on
% for i = 1:length(Uinf_)
%     plot(real(p_values_collect(i,:)),imag(p_values_collect(i,:)))
% end
% hold off

figure()
plot(Uinf_,p_values)
grid on
grid minor
ylabel("$max(Re(p_i))$",'Interpreter','latex')
xlabel("$U_{\infty}$",'Interpreter','latex')
hold off

end





