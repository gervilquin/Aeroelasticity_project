%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aeroelastic model of a Pazy wing
% Authors: Gerard Villalta & Antoni Alberto
% Date: 18/11/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc, clear, close all
addpath(genpath('functions'))

% Case solution active/deactivate
solve_static = true; 
solve_diverge = true;
solve_modal = true;
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
U_inf = 30;  % Freestream velocity (m/s)
AoA = 5;     % Wing angle of attack (º)
rho_inf = 1.3; % Reference air density (kg/m^3)

%% 2. Mesh construction
% Mesh definition
y_el = ComputeYcoordinates(y_sec,Nsec,Nesec);  % Elements' nodes Y coordinate (m)
Nnodes = length(y_el);
Nel = Nnodes - 1;

% Connectivity matrices definition
Tn = ConnectivityElements(Nel);
Ts = ConnectivitySubsets(Tn,Nsec,Nesec);

% Dirichlet boundary conditions
Up = [  deg2rad(AoA) 1 1;
        0 1 2;
        0 1 3];

%% 3. Structural modelling
% Define stiffness matrix
K = ComputeKmatrix(y_el,Tn,EI,GJ);
% Define mass matrix
M = ComputeMmatrix(y_el,Tn,Ts,material,t,c,x_sc*c,RibExist);

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
[I_au_0, I_au_1, I_au_2] = ComputeDisplacementsCoupling(y_el,Tn,U_inf);
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
    saveas(gcf,'report/figures/static_solution','epsc')
end


%% 7. Static solution with aerocoupling
if solve_static == true
    [Meff,Ceff,Keff] = ComputeEffectiveMatrix(K,M,y_el,Tn,U_inf,rho_inf,A,S,I_fL);
    AeroStaticSolver(y_el,u_static,If,Ip,I_fL,I_au_0,S,A,Keff,U_inf,rho_inf,AoA,true);
    saveas(gcf,'report/figures/Aeroelastic_solution','epsc')
end

%% 8. Compute divergence
if solve_diverge == true
    Uinf_ = linspace(0.1,150,50);
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
    ylabel("$w_{sc}^{tip}$","Interpreter","latex")
    grid on
    grid minor
    legend([strcat("Divergence $U_{\infty}$ = ",string(round(U_diverg(end),2))," m/s")],'location','northwest',"Interpreter","latex")
    hold off
    saveas(gcf,'report/figures/divergence_velocity','epsc')
end

%% 9. Modal analysis
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
    ylabel("$\theta$","Interpreter","latex")
    
    
    %Plot mode shapes in the vertical displacement 
    nexttile(tcl)
    hold on
    for i =1:Nm
        plot(y_el,Phi(2:3:end,i))
    end
    grid on
    grid minor
    hold off
    ylabel("$w_{sc}$","Interpreter","latex")

    %Plot the mode shapes in deflection
    nexttile(tcl)
    hold on
    for i =1:Nm
        plot(y_el,Phi(3:3:end,i),'DisplayName',strcat("Mode ",string(i),": ",string(round(freq(i),0))," Hz"))
    end
    grid on
    grid minor
    hold off
    ylabel("$\gamma$","Interpreter","latex")
    xlabel("Y [mm]",'Interpreter','latex')

    % Add legend to the side
    Lgnd = legend('show','interpreter','latex');
    Lgnd.Layout.Tile = 'East';

    saveas(gcf,'report/figures/Natural_modes','epsc')

end



%% 10. Flutter solver
if solve_flutter == true
%Uinf_ = logspace(-10,-1,100);
%Uinf_ = linspace(0.1,U_diverg(end),3);
Uinf_ = linspace(0.1,60,100);
%Uinf_ = [20];

% Get the eigenvalues of M and K
N_reduced_ = [3];

for k = 1:length(N_reduced_)

N_reduced = N_reduced_(k);
%N_reduced = 3;
[Vr, Dr] = eigs(K(If,If),M(If,If),N_reduced,'sm');

% Initialize matrices
p_values = zeros(length(Uinf_),1);
p_values_collect = zeros(length(Uinf_),2*N_reduced);

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
        -1*eye(size(Keff_red,1)) zeros(size(Keff_red,1))];

    % Compute eigen values
    [Vd, Dd] = eigs(D,eye(length(D)),length(D),'lm');

    p_values(i) = max(real(-1./diag(Dd)));
    p_values_collect(i,:) = -1./diag(Dd);
end


%% Plots

figure()
hold on
for i = 1:length(p_values_collect(1,:))/2
    mode = (i-1)*2;
    plot(real(p_values_collect(:,mode+1)),imag(p_values_collect(:,mode+1)),'DisplayName',strcat("Mode ",string(i)," positive"))
    plot(real(p_values_collect(:,mode+2)),imag(p_values_collect(:,mode+2)),'DisplayName',strcat("Mode ",string(i)," negative"))
    plot(real(p_values_collect(1,mode+1)),imag(p_values_collect(1,mode+1)),'HandleVisibility','off','Marker','o','Color','red')
    plot(real(p_values_collect(1,mode+2)),imag(p_values_collect(1,mode+2)),'HandleVisibility','off','Marker','o','Color','red')
    plot(real(p_values_collect(end,mode+1)),imag(p_values_collect(end,mode+1)),'HandleVisibility','off','Marker','x','Color','g')
    plot(real(p_values_collect(end,mode+2)),imag(p_values_collect(end,mode+2)),'HandleVisibility','off','Marker','x','Color','g')
end
xline(0,'color','k','HandleVisibility','off')
yline(0,'color','k','HandleVisibility','off')
grid on
grid minor
xlabel("Re($\hat{p}_i$)",'Interpreter','latex')
ylabel("Im($\hat{p}_i$)",'Interpreter','latex')
legend("Location",'eastoutside','Interpreter','latex')
saveas(gcf,'report/figures/p_value_complex','epsc')

hold off

figure()

plot(Uinf_,p_values)
grid on
grid minor
ylabel("$max{Re(p_i)}$",'Interpreter','latex')
xlabel("$U_{\infty}$",'Interpreter','latex')
yline(0,'color','k')
legend([strcat("$U_{flutter}$ = ",string(round(compute_flutter_velocity(p_values,Uinf_),2))," m/s")],'Interpreter','latex')
hold off
saveas(gcf,'report/figures/Flutter_velocity','epsc')

end

end





