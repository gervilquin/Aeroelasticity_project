clear
close all;

% Data
a_ = -0.2
r_ = 0.5
x_ = 0.3
mu_ = 20
sig_ = 0.01:001:1.8;
beta = 0.5; %<- for newton raphson

% Matrices 
M = @(mu,r,x) mu*[r^2. -x; -x, 1]; %> fast way to define functions
K =@(mu,r,sig) mu
M1 = @(a) [1/8+a^2, a; a,1];
C1
K1
C2

% Tjerodorzem function 

%       x^4     x^3     x^2     x       1
pr = [  0.5,    0,      0.0765, 0       1.8632e-4];
p1 = 
q 
F = @(k) polyval(pr,k)/polyval(q,k)
G 
C 
dFdk
dGdk 
dCdk

% Loop through

for i = 1:length(sig_)
    
    % evaluate the matrices that have been defined
    M_ = M(mu_, r_, x_);
    K_ = k(mu_,r_,sig_(1));
    M1_ = M1(a_);
    C1_ = C1(a_);
    K1_ = K1(a_);
    C2_ = C2(a_);

    % system matric
    D = @(k) K_\(M_+M1_+1i/k*(C2_ - C(k)*C1_)+1/k^2*C(k)*K1_);
    dDdk = @(k) K_\((dFdk(k)/k^2 - 2*F(k)/k^3)*K1_ + (dGdk(k)/k - G(k)/k^2)*C1_) ...
        +1i*K_\((dGdk(k)/k^2 - 2*G(k)/k^3)*K1_ + (dFdk(k)/k - F(k)/k^2)*C1_ - 1/k^2*C2_);

    %% METHOD 2 (inital Guess)
    % Trial values f the reduced frequeccny

    k_test = linspace(0.1,1,10);
    %k_test = logspace(-3,3,500);% by increasing the number of points this will lead to a better solution with method 2
    for j = 1:length(k_test)
        % eValuate the system matrix
        D_ = D(k_test(j));

        %det(D(k) - lamnda*I) = 0 ->  eigenvalues problem
        %[~,L ] = eig(D_); % if this function is used to compute the eigen
        %values, then this may return the egenvalues unsorted and will lead
        %to jumps in the solution, it is required to verify that this is
        %not an issue
        c = [1,-(D_(1,1) + D_(2,2)), D_(1,1)*D_(2,2)-D_(1,2)*D_(2,1)];

        % find the roost of the 2nd order polynomial

        l_test(1) = -c(2)/(2*c(1)) + 1/(2*c(1))*sqrt(c(2)^2 - 4*c(1)*c(3));
        l_test(1) = -c(2)/(2*c(1)) - 1/(2*c(1))*sqrt(c(2)^2 - 4*c(1)*c(3));

        % Get artifical damping (g) and freestreem velocity (U)
        for k=1:2
            g(i,k) = imag(l_test(k))/real(l_test(k));
            U(j,k) = 1/(k_test(j)*sqrt(real(l_test(k))));
        end

    end

    % look for sign changes

    kF = [];
    UF = [];

    for k = 1:2
        j0 = find(g(1:end-1,k).*g(2:end,k) < 0);

        %if ~isempty(j0) % adding this if is redundant, the for loop does
        %it alone
            for j = 1:length(j0) % ignored if the j0 vector is empty
                kF = [kF,interp1(g([j0(j),j0(j) +1],k), k_test([j0(j),j0(j)+1]),0)]; % implement the linear interpolation between the two extrems of the sign change
                UF = [UF,interp1(g([j0(j),j0(j) +1],k), U([j0(j),j0(j)+1],k),0)];
            end
        %end
    end
    Uf(i) = UF(UF == min(UF));
    kf(i) = kF(UF == min(UF));

    % plot the results  as shon in the slides
    % >> plot code lines <<

    %%  METHOD 1: Newton Raphson
    
    %First gess
    x = [(1/(Uf(i) * kf(i)))^2;kf(i)];
    tol = sum(x.^2);
    j = 0;

    % loop 
    while(abs(tol) > 10^-8 &&  j < 500
        % update the iteration
        j = j+1;
        %values of lamnda and k
        % evalaute ystem matrix and tis derivate
        % evaluate determinat component (real and imaginary)

        % update solution
        dx = -J\f;
        x = x +beta*dx;
        % compute convergence
        tol = sum(f.^2);

    end
    Uf_(i)= 1/(x(2)*sqrt(x(1)));
    
end

