clear
close all
clc
fprintf('Welcome to Code for Aeroelasticity Task\n')
fprintf('CAT v20230107\n')
fprintf('Analysis:  Dynamic')
fprintf("\n |\\__/,|   (`\\\n |_ _  |.--.) )\n ( T   )     /\n(((^_(((/(((_/\n")

% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));


%% 1 - INPUT DATA
mode = "flutter";
u_inf=linspace(1,60,100);  % freestream velosidá
% u_inf = 10;
alpha0 = deg2rad(0);
n_reduction      = 12;        % number of modes for the order reduction

% 1.1 - Structural data
try load("stiffness/estimated.mat")
catch, error("Please create stiffness/estimated.mat by running main_stiffness first!") 
end

xsc             = estimated.x_sc;     % x coordinates of the shear center
GJ              = estimated.GJ;       % torsional stifness
EI              = estimated.EI;       % bending stifness
rho_al          = 2795;       % aluminium density
rho_po          = 930;        % polymer density
nodes_subset    = 12;                  % nodes per subseb
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
if strcmp(mode,"flutter")
    fprintf('u_inf(min): %f [m/s]\n',min(u_inf))
    fprintf('u_inf(max): %f [m/s]\n',max(u_inf))
else
    fprintf('u_inf:      %f [m/s]\n',u_inf)  
end
fprintf('n_reduction:%d \n',n_reduction)


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
    %        y      -> global position of the nodes y_ribs -> global
    %        polition of the ribs l      -> element length
    %     y(i)----------y(i+1)----------y(i+2)...
    %             l(e)          l(e+1)
    [y, y_ribs, l] = compute_nodes(ns);
    % 2.1.5 - Assembly of global matrices
    %        DOFs   -> Total number of DOFs K      -> bar{K} in the slides,
    %        global stifness matrix M      -> bar{M} in the slides, global
    %        mass matrix
    [DOFs, Kbar, Mbar] = Assembly(y, l,GJ,EI,rho_al,rho_po,Tn2,ns,xsc,Ts,myNACA,hr);   
    if n_reduction<=0, n_reduction = DOFs+n_reduction; end
    check_simmetry(Kbar,'K structural');
    check_simmetry(Mbar,'M structural');
    % 2.1.6 - Solve eigenproblem (from comutational engineering), idk if it
    %         needs to be solved for n=DOFs-prescribed_dofs
    [phi,w2]  = compute_modal(Kbar,Mbar,n_reduction,DOFs,4:DOFs);
    if strcmp(mode,"modal"), plot_modal(y,phi,w2,n_reduction); end
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
    [collocation_points, aerodynamic_point,horse_shoe] = ....
        aero_disc(y,c,xsc,0);
    % 2.2.2 - Obtain the normals to the surfaces
    [Sf, n] = compute_normals_surface(horse_shoe,c); 
    % 2.2.3 - Obtain AIC matrix
    A = compute_A(collocation_points,horse_shoe,n);
    % 2.2.4 - Obtain V matrix
    %         by comparison of slide 22 project with slide 1 S06 I
    %         understand that V has this value [A]{Gamma} = -Uinf{alpha}
    %         <-> [A]{Gamma} = [V(Uinf)]{alpha}
    V = @(U_inf) -eye(length(aerodynamic_point),length(aerodynamic_point)).*U_inf;
    % 2.2.5 - Obtain S matrix
    %         by comparison of slide 22 project with slide 1 S06 I
    %         understand that S has this value L =
    %         rhoinf*Uinf*l[i]*Gamma[i] <-> {L} = [S(Uinf)]{Gamma} be
    %         awarned that S is a func of u_inf as proposed in S06 - slide
    %         10
    S = compute_S(rho_inf,l);

% 2.3 - COUPLING MODULE
    % 2.3.1 - Obtain I_alpha_u using shape functionsç
    %         be awarned that I_alpha_u0, I_alpha_u1 and I_alpha_u2 are in
    %         general funcs of U_inf, as proposed in S06 - slide 6 and
    %         discussion topic "Definition of I_alpha_u matrices" in the
    %         forum of the course
    [I_alpha_u0, I_alpha_u1,I_alpha_u2] = struct2aero(y,aerodynamic_point,DOFs,xsc,xac);
    % 2.3.1 - Obtain I_f_L dividing lift between neighbour structural nodes
     %        maybe the error comes from here?
    [I_f_L]     = aero2struct(y,aerodynamic_point,DOFs,xsc,xac);


if strcmp(mode,"flutter")
    %% 3 - LOOP OVER VELOCITIES
    for n = 1:length(u_inf)
        % 3.1 - EVALUATE AND COMPUTE NEEDED MATRICES
            % 3.1.2 - Evaluate aerodynamic matrices and I_alpha_u
            %         as proposed in S06 - slide 15
            S_now = S(u_inf(n));
            V_now = V(u_inf(n));
            I_alpha_u0_now = I_alpha_u0(u_inf(n));
            I_alpha_u1_now = I_alpha_u1(u_inf(n));
            I_alpha_u2_now = I_alpha_u2(u_inf(n));
            % 3.1.2 - Compute system matrices
            %         as proposed in S06 - slide 10
            M_a = -I_f_L*S_now*inv(A)*V_now*I_alpha_u2_now;
            C_a = -I_f_L*S_now*inv(A)*V_now*I_alpha_u1_now;
            K_a = -I_f_L*S_now*inv(A)*V_now*I_alpha_u0_now;
            M_eff = Mbar + M_a;
            K_eff = Kbar + K_a;
            C_eff = C_a;
            % 3.1.3 - Apply order reduction
            %         as proposed in S06 - slide 12
%             [phi,w2]  = compute_modal(K_eff,M_eff,n_reduction,DOFs,4:DOFs);
            mi = 1;
            M_eff_tilde  = phi(:,mi:end)'*M_eff*phi(:,mi:end);
            C_eff_tilde  = phi(:,mi:end)'*C_eff*phi(:,mi:end);
            K_eff_tilde  = phi(:,mi:end)'*K_eff*phi(:,mi:end);
        % 3.2 - SOLVE THE SYSTEM
            % 3.2.1 - Compute D_prime
            %         as proposed in S06 - slide 13, method b1)
            D_prime = compute_D_prime(K_eff_tilde,C_eff_tilde,M_eff_tilde);
            % 3.2.2 - Solve Eigenproblem
            lambda_p = solve_system(D_prime);
            % 3.2.3 - Find maximum real part of lambda
            p(n)          = max(real(-1./lambda_p)); % min¿?
            p_real(n,:)   = real(-1./lambda_p); 
            p_im(n,:)     = imag(-1./lambda_p);

%              figure 
%             hold on
%             grid minor
%             plot(p_real(n,:) , p_im(n,:));
%             xlabel('p_real')
%             ylabel('p_im')


    end
    
     
    %% 4 - RESULTS
    
    % 4.1 - Interpolate U_F
    %       for g = 0;
    %         ???????
    % 4.2 - Plot g - u_inf
    figure 
    hold on
    grid minor
    plot(u_inf,p,'LineWidth',1.5)
    xlabel('u_\infty')
    ylabel('p')

%     figure 
%     hold on
%     grid minor
%     plot(p_real,p_im);
%     xlabel('u_\infty')
%     ylabel('p')

%     figure 
%     hold on
%     grid minor
% %     plot(u_inf,p_real);
%     xlabel('u_\infty')
%     ylabel('p')

end

%% APPENDIX A - PLOTS
% plot_lift(L,collocation_points);
% draw_model(y,y_ribs,myNACA,alpha,Tn,xsc,delta,4e-3,collocation_points,
% aerodynamic_point,horse_shoe,gamma);
% plot_struct_disc(y,y_ribs,myNACA,alpha,Tn,xsc,delta,4e-3);


%% APPENDIX B -DORAEMON 
% L = compute_L(rho_inf,u_inf,S,alpha,A); alpha =
% ones(size(y)).*deg2rad(10); % Hardcode, just to see beautiful plots delta
% = ones(size(y)).*1e-4;             % Hardcode, just to see beautiful
% plots gamma = zeros(size(delta)); theta = ones(size(delta)).*0.3;
 % 2.1.10- Mass-normalization of the eigenmodes (slide 12 S06)
%     phi      = mass_normalization(phi,M);

%% APPENDIX C - POSSIBLE FAILURES
% 1. solver is not implemented? this has been so easy 
% 2. does C exist? 
% 3.under matrix of the prime is not an eye 
% 4. Stillproblems with the last secion? 
% 5. We de not reduce but there are still problems 
% 6. Nodes_subset changes the order of magnitude of the results?
% 7. Is I_f_L well calculated? We think so, f is coherent with L 
% 8.Modal analysis does not result: causes?
    % 8.1 - Mass matrix is wrong? OR in static is well
    % 8.2 - order_reduction is wrong, OR in static is well
%
