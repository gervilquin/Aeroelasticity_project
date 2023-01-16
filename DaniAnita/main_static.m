clear
clc
close all

%% 1 - INPUT DATA
fprintf('Welcome to Code for Aeroelasticity Task\n')
fprintf('CAT v20230107\n')
fprintf('Analysis:   static')
fprintf("\n |\\__/,|   (`\\\n |_ _  |.--.) )\n ( T   )     /\n(((^_(((/(((_/\n")
tic

% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

% 1.0 - Mode (normal, divergence, debug)
mode = 'divergence';
u_inf=linspace(1,180,100);  % freestream velosidá
% u_inf = 55;  % freestream velosidá
n_reduction     = 5;       % number of modes for the order reduction
alpha0 = deg2rad(7);

% u_inf = 50;
% 1.1 - Structural data
try load("stiffness/estimated.mat")
catch, error("Please create stiffness/estimated.mat by running main_stiffness first!")
end
xsc             = estimated.x_sc;     % x coordinates of the shear center
GJ              = estimated.GJ;       % torsional stifness
EI              = estimated.EI;       % bending stifness
rho_al          = 2795;       % aluminium density
rho_po          = 930;        % polymer density
nodes_subset    = 12;         % nodes per subset
nodes_rib       = 18;         % nodes in rib
dead_ribs       = 3;          % ribs with no nodes between them and the next rib
c               = 100e-3;     % chord
hr              = 4e-3;       % thickness of the ribs, m
prescribed_dofs = 3;          % prescribed dofs, must be the first ones 
% 1.2 - Aerodynamics data

xac = c/4;
rho_inf = 1.23;           % freestream flow density
% Wing profile
t = 0.018; % Harcode, this must be a 0.18 but then the profile is fat
myNACA = @(x) 5.*t.*(0.2969.*sqrt(x/c)-0.1260.*(x/c)-0.3516.*(x/c).^2+...
       0.2843.*(x/c).^3-0.1015.*(x/c).^4);
fprintf('Mode:       %s\n',mode)
fprintf('AOA:        %f [deg]\n',rad2deg(alpha0))
if strcmp(mode,'divergence')
    fprintf('u_inf(min): %f [m/s]\n',min(u_inf))
    fprintf('u_inf(max): %f [m/s]\n',max(u_inf))
else
    fprintf('Velocity:   %f [m/s]\n',u_inf)
end
%% 2 - PREPROCESS/INITIALIZATION

% 2.1 - STRUCTURAL MODULE
    % 2.1.1 - Create ns 
    %        ns -> number of local nodes per element
    ns = compute_ns(nodes_subset,nodes_rib,dead_ribs);
    % 2.1.2 - Create T 
    %        Ts -> nodal connectivities matix
    Ts = compute_Ts(ns);
    % 2.1.3 - Create Tn
    %        Tn -> element connectivities matrix
    [Tn,Tn2] = compute_Tn(ns);
    % 2.1.4 - Create structural mesh
    %        y      -> global position of the nodes
    %        y_ribs -> global polition of the ribs
    %        l      -> element length
    %     y(i)----------y(i+1)----------y(i+2)...
    %             l(e)          l(e+1)
    [y, y_ribs, l] = compute_nodes(ns);
    % 2.1.5 - Assembly of global matrices
    %        DOFs   -> Total number of DOFs
    %        K      -> bar{K} in the slides, global stifness matrix
    %        M      -> bar{M} in the slides, global mass matrix
    [DOFs, Kbar, Mbar] = Assembly(y, l,GJ,EI,rho_al,rho_po,Tn2,ns,xsc,Ts,myNACA,hr);
    % 2.1.6 - Solve eigenproblem (from comutational engineering)
    [phi,w2]  = compute_modal(Kbar,Mbar,n_reduction,DOFs,4:DOFs);
    % 2.1.7- Create vector omega2 as in slide 14 S06
    Omega2    = diag(w2);
% 2.2 - AERODYNAMIC MODULE   
    % 2.2.1 - Create aerodynamics mesh
    %        collocation_points -> Points marked with a x in the slides
    %        aerodynamic_point  -> Points marked with a o in the slides
    %        horse_shoe         -> Struct containing matrices A,B,C,D each
    %                              one containing the position of one of
    %                              the points of each horseshoe
    % z0 y->   B--o--C
    %  x       |     |
    %  |       |  x  |
    %  v       |     |
    %          |     |
    %          |     |
    %          A-----D
    [collocation_points, aerodynamic_point,horse_shoe] = ...
        aero_disc(y,c,xsc,0);
    % 2.2.2 - Obtain the normals to the surfaces
    [Sf, n] = compute_normals_surface(horse_shoe,c); 
    % 2.2.3 - Obtain AIC matrix
    A = compute_A(collocation_points,horse_shoe,n);
    % 2.2.4 - Obtain V matrix
    %         by comparison of slide 20 project with slide 1 S06 
    %         I understand that V has this value 
    %         [A]{Gamma} = -Uinf{alpha} <-> [A]{Gamma} = [V(Uinf)]{alpha}
    V = @(U_inf) -eye(length(aerodynamic_point),length(aerodynamic_point)).*U_inf;
    % 2.2.5 - Obtain S matrix
    %         by comparison of slide 22 project with slide 1 S06 
    %         I understand that S has this value 
    %         L = rhoinf*Uinf*l[i]*Gamma[i] <-> {L} = [S(Uinf)]{Gamma}
    %         be awarned that S is a func of u_inf as proposed in S06 -
    %         slide 10
    S = compute_S(rho_inf,l);

% 2.3 - COUPLING MODULE
    % 2.3.1 - Obtain I_alpha_u using shape functionsç
    %         be awarned that I_alpha_u0, I_alpha_u1 and I_alpha_u2 are 
    %         in general funcs of U_inf, as proposed in S06 - slide 6 and 
    %         discussion topic "Definition of I_alpha_u matrices" in the 
    %         forum of the course
    [I_alpha_u0, I_alpha_u1,I_alpha_u2] = struct2aero(y,aerodynamic_point,DOFs,xsc,xac);
    % 2.3.1 - Obtain I_f_L dividing lift between neighbour structural nodes
     %        maybe the error comes from here?
    [I_f_L]     = aero2struct(y,aerodynamic_point,DOFs,xsc,xac);

%% 3 - STATIC CASE


if strcmp(mode,'debug')
     dbstop in main_static at 160
end

n         = 1; % to loop over the velocities        
breakflag = 0; % if we need to break the loops

while 1
    % 3.1 - INITIALIZE ALPHA AND U
    alpha = zeros(length(aerodynamic_point),1);
    u = zeros(DOFs,1);
    % 3.2 - COMPUTE OUT-LOOP MATRICES
    S_ = S(u_inf(n));
    V_ = V(u_inf(n));
    I_alpha_u0_ = I_alpha_u0(u_inf(n));
    I_alpha_u1_ = I_alpha_u1(u_inf(n));
    I_alpha_u2_ = I_alpha_u2(u_inf(n));
    % 3.3 - CONVERGE THE SYSTEM
    nite = 0; % number of iterations
    while 1
        nite = nite +1;
        % 3.3.1 - Compute the lift force
        L = S_*inv(A)*V_*...
            (alpha+ones(length(aerodynamic_point),1).*alpha0);        
        % 3.3.2 - Obtain cl
        CL(n) = sum(L)/(0.5*rho_inf*u_inf(n)^2*sum(Sf));
        cl = L./(0.5*rho_inf*u_inf(n)^2*sum(Sf)); 
        % 3.3.3 - Convert lift to structural forces
        f = I_f_L*L;        
        % 3.3.4 - Solve the system
        u_np1 = reduced_order(Kbar,Mbar,n_reduction,DOFs,4:DOFs,f);
        u_np1 = (u_np1);
        % 3.3.5 - Obtain the INDUCED alpha by the struct deformations
        alpha_np1 = I_alpha_u0_*u_np1;           
        % 3.3.6 - If in debug mode, check lift and open vars panel
        if strcmp(mode,'debug')
            L2 = compute_L(rho_inf,u_inf(n),l,...
                    alpha+ones(length(aerodynamic_point),1).*alpha0,A);
            if max(abs(L-L2))>1e-2, error("Different L!"); end
            plot_vars_static(aerodynamic_point,y,u_np1,f,alpha_np1,L);
            
         end                  
        % 3.3.7 - Check if the problem has diverged by the nite
        if nite>10e2
            fprintf("CAT diverged for u_inf=%f\n",u_inf(n))
            fprintf("end algorithm and print last ite\n")
            breakflag = 1;
            break;
        end 
        % 3.3.8 - Check if the problem has converged in the struct
        if max(max(abs(u-u_np1)))<1e-6, break; end
        % 3.3.9 - Update values
        u     = u_np1;
        alpha = alpha_np1;         
    end

    % 3.4 - SAVE THE DESIRED RESULTS
        % 3.4.1 - Save the deformations at tip
        w_tip(n)      = u(end-1);
        theta_tip(n)  = u(end-2);
        gamma_tip(n)  = u(end);
    
    % 3.5 - END?
    %       If we have reached the end or diverged, exit this loop
    if n==length(u_inf) || breakflag, break; end
    n = n + 1;
    
end
et = toc;
fprintf("--CAT RESULTS--\n")
fprintf("Number of iterations = %d\n",nite);

if strcmp(mode,'normal')
        plot_deformations(y,u); 
        sca   = 1;
        alpha = sca*u(1:3:end);% + ones(length(y),1).*alpha0;
        delta = sca*u(2:3:end);
        gamma = sca*u(3:3:end);
        draw_model(y,y_ribs,myNACA,alpha,Tn,xsc,delta,hr,...
                    collocation_points, aerodynamic_point,horse_shoe,gamma);
elseif strcmp(mode,'divergence')
    figure
    hold on
    grid minor
    plot(u_inf(1:n),w_tip,'LineWidth',1.5)
    xlabel('u_\infty')
    ylabel('w_{tip}')
    
    figure
    hold on
    grid minor
    plot(u_inf(1:n),rad2deg(theta_tip),'LineWidth',1.5)
    xlabel('u_\infty')
    ylabel('\Theta_{tip}')

    figure
    hold on
    grid minor
    plot(u_inf(1:n),rad2deg(gamma_tip),'LineWidth',1.5)
    xlabel('u_\infty')
    ylabel('\gamma_{tip}')
    fprintf("Divergence mode, printing max values\n")
end
fprintf("w at wing tip        = %f [mm]\n",w_tip(end)*1e3);
fprintf("theta at wing tip    = %f [deg]\n",rad2deg(theta_tip(end)));
fprintf("gamma at wing tip    = %f [deg]\n",rad2deg(gamma_tip(end)));
fprintf("CL calculated        = %f\n",CL(end));
fprintf("CL diff              = %f\n",CL(end)-2*pi*alpha0)
fprintf("Elapsed time         = %f [s]\n",et);
if strcmp(mode,"divergence")
    fprintf("Divergence mode, printing mean values\n")
    fprintf("w at wing tip        = %f [mm]\n",mean(w_tip)*1e3);
    fprintf("theta at wing tip    = %f [deg]\n",mean(rad2deg(theta_tip)));
    fprintf("gamma at wing tip    = %f [deg]\n",mean(rad2deg(gamma_tip)));
    fprintf("CL calculated        = %f\n",mean(CL));
    fprintf("CL diff              = %f\n",mean(CL)-2*pi*alpha0)
    fprintf("Elapsed time         = %f [s]\n",et);
end
%% APPENDIX A - PLOTS

% plot_lift(L,collocation_points);
% draw_model(y,y_ribs,myNACA,alpha,Tn,xsc,delta,4e-3,collocation_points, aerodynamic_point,horse_shoe,gamma);
% plot_struct_disc(y,y_ribs,myNACA,alpha,Tn,xsc,delta,4e-3);


%% APPENDIX B -DORAEMON 
% L = compute_L(rho_inf,u_inf,S,alpha,A);
% alpha = ones(size(y)).*deg2rad(10); % Hardcode, just to see beautiful plots
% delta = ones(size(y)).*1e-4;             % Hardcode, just to see beautiful plots
% gamma = zeros(size(delta));
% theta = ones(size(delta)).*0.3;
 % 2.1.10- Mass-normalization of the eigenmodes (slide 12 S06)
%     phi      = mass_normalization(phi,M);

%% APPENDIX C - POSSIBLE FAILURES
% 1. solver is not implemented? this has been so easy
% 2. does C exist?
% 3. I_f_L does not give much confidence. But it is checked
% 4. under matrix of the prime is not an eye
% 5. Why the reducing changes the results so much??? suspicious
% 6. We can try the static case. Why is it logical but the results for the
%    dynamic case are not?
% 7. Still problems with the last secion?
% 8. We de not reduce but there are still problems
% 9. Nodes_subset changes the order of magnitude of the results?
% 10. Does the problem comes from the assembly?





