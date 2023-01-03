clear
close all

% Data
a_ = -0.2; %shear center location
r_ = 0.5; 
x_ = 0.3;
mu_ = 20; % Mass/density ratio
sig_ = 0.01:0.01:1.8; % Stiffness ratio
beta = 0.5;

% Matrices definitions
M = @(mu,r,x)   mu*[r^2, -x; -x,     1];
K = @(mu,r,sig) mu*[r^2,  0;  0, sig^2];
M1 = @(a) [1/8+a^2, a; a, 1];
C1 = @(a) [2*a^2-1/2, 2*a+1; 2*a-1, 2];
K1 = @(a) [1+2*a, 0; 2, 0];
C2 = @(a) [a-1/2, 0; 1, 0];

% Theodorsen's function
%        x^4     x^3     x^2           x          1
pr = [   0.5,      0, 0.0765,          0, 1.8632e-4];
pi = [        -0.108,      0, -8.8374e-4,         0];
q  = [     1,      0, 0.0921,          0, 1.8632e-4];
F = @(k) polyval(pr,k)/polyval(q,k);
G = @(k) polyval(pi,k)/polyval(q,k);
C = @(k) F(k) + 1i*G(k);
dFdk = @(k) polyval(polyder(pr),k)/polyval(q,k) - polyval(polyder(q),k)/polyval(q,k)*F(k);
dGdk = @(k) polyval(polyder(pi),k)/polyval(q,k) - polyval(polyder(q),k)/polyval(q,k)*G(k);
dCdk = @(k) dFdk(k) + 1i*dGdk(k);

% Loop through parameter
for i = 1:length(sig_)
    
    % Evaluate
    M_ = M(mu_,r_,x_);
    K_ = K(mu_,r_,sig_(i));
    M1_ = M1(a_);
    C1_ = C1(a_);
    K1_ = K1(a_);
    C2_ = C2(a_);
    
    % System matrix
    D = @(k) K_\(M_+M1_+1i/k*(C2_-C(k)*C1_)+1/k^2*C(k)*K1_);
    
    dDdk = @(k) K_\((dFdk(k)/k^2 - 2*F(k)/k^3)*K1_ + (dGdk(k)/k - G(k)/k^2)*C1_) ...
           + 1i*K_\((dGdk(k)/k^2 - 2*G(k)/k^3)*K1_ + (dFdk(k)/k - F(k)/k^2)*C1_ - 1/k^2*C2_);
    
    %% Method 2 (initial guess)
    % Trial values of the reduced frequency
    k_test = linspace(0.1,1,10);
    
     for j = 1:length(k_test)
        % Evaluate the system amtrix
        D_ = D(k_test(j));
        % det(D(k) - lambda*I) = 0
        %[~,L] = eig(D_);
        c = [ 1, -(D_(1,1)+D_(2,2)), D_(1,1)*D_(2,2)-D_(1,2)*D_(2,1)];
        % Find the roots of the 2nd order polynomial (=eigenvalues)
        l_test(1) = -c(2)/(2*c(1)) + 1/(2*c(1))*sqrt(c(2)^2-4*c(1)*c(3));
        l_test(2) = -c(2)/(2*c(1)) - 1/(2*c(1))*sqrt(c(2)^2-4*c(1)*c(3));
        % get artificial damping (g) and freestream velocity (U)
        for k=1:2
            g(j,k) = imag(l_test(k)/real(l_test(k)));
            U(j,k) = 1/(k_test(j)*sqrt(real(l_test(k))));
        end
    end
    % Look for sign changes
    kF = [];
    UF = [];
    for k = 1:2
        j0 = find(g(1:end-1,k).*g(2:end,k)<0); % Multiply each result, if they have different sign, the result will be negative
        for j = 1:length(j0)
            kF = [kF,interp1(g([j0(j),j0(j)+1],k),k_test([j0(j),j0(j)+1]),0)];
            UF = [UF,interp1(g([j0(j),j0(j)+1],k),U([j0(j),j0(j)+1],k),0)];
        end
    end
    Uf(i) = UF(UF==min(UF));
    kf(i) = kF(UF==min(UF));
    
    % Plot stability curves
    if sig_(i)==1
        figure(1)
        yyaxis left
        cla
        hold on
        grid on
        plot(k_test,g(:,1));
        plot(k_test,g(:,2));
        plot(kF,zeros(size(kF)),'xb','markersize',10);
        ylabel('g');
        yyaxis right
        cla
        hold on
        plot(k_test,U(:,1));
        plot(k_test,U(:,2));
        plot(kF,UF,'xr','markersize',10);
        ylim([0,4]);
        ylabel('U_{F}/\omega_{\theta}b')
        set(gca,'xscale','log');
        xlabel('\kappa_F')   
        drawnow
    end
    
    %% Method 1 (Newton-Raphson) 
    
%     x = [];
%     while abs(tol) > 1e-8 && j < 500
%         % Update iteration
%         j = j+1;
%         % Values of lambda and k
%         l = x(1);
%         k = x(2);
%         % Evaluate system matrix and its derivative
%         D_ = D(k);
%         Dr = real(D_);
%         Di = iamg(D_);
%         dDdk_ = dDdk(k);
%         Drk = real(dDdk_);
%         DiK = imag(dDdk_);
%         % Evalute determinant components (real and imaginary)
%         
%         % Evaluate jacobina matrix
%         J 
%         % Update solution
%         dx = -J/f;
%         x = x + beta*dx;
%         % Compute convergence
%         
% 
%     end
    % Sotee flutter speed
    %Uf_(i) = 1/(x(2)*sqrt(x(1)));
    
end

figure(2)
cla
hold on
box on
grid on
plot(sig_,Uf,'ok')
%plot(sig_,Uf_,'-b')
xlim([0,1.8]);
ylim([1.5,4]);
xlabel('\sigma');
ylabel('U_F/\omega_{\theta}b');
