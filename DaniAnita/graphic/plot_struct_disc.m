function plot_struct_disc(y,y_ribs,myNACA,alpha,Tn,xsc,delta,thickness_rib)
    
    % Density quantity of points for drawing the NACAs
    N=30;
    
    % X position for points for NACAS
    x_for_NACAs = linspace(0,0.1,N);
    
    figure
    hold on
    
    % Plot the structural nodes and the shear center line
    plot3(xsc.*ones(size(y)),y,delta,'rs-');
    for i=1:length(y)
        A(i,:) = [0,y(i),0];
        B(i,:) = [7.5e-3,...
                  y(i),...
                  myNACA(7.5e-3)+delta(i)];
        C(i,:) = [7.5e-3,...
                  y(i),...
                  -myNACA(7.5e-3)+delta(i)];
        D(i,:) = [35e-3,y(i),2.5e-3+delta(i)];
        E(i,:) = [35e-3,y(i),-2.5e-3+delta(i)];
        F(i,:) = [65e-3,y(i),2.5e-3+delta(i)];
        G(i,:) = [65e-3,y(i),-2.5e-3+delta(i)];
        H(i,:) = [95e-3,...
                  y(i),...
                  myNACA(95e-3)+delta(i)];
        I(i,:) = [95e-3,...
                  y(i),...
                  -myNACA(95e-3)+delta(i)];
        J(i,:) = [100e-3,...
                  y(i),...
                  delta(i)];
        A(i,:) = A(i,:) - [xsc,y(i),delta(i)];
        B(i,:) = B(i,:) - [xsc,y(i),delta(i)];
        C(i,:) = C(i,:) - [xsc,y(i),delta(i)];
        D(i,:) = D(i,:) - [xsc,y(i),delta(i)];
        E(i,:) = E(i,:) - [xsc,y(i),delta(i)];
        F(i,:) = F(i,:) - [xsc,y(i),delta(i)];
        G(i,:) = G(i,:) - [xsc,y(i),delta(i)];
        H(i,:) = H(i,:) - [xsc,y(i),delta(i)];
        I(i,:) = I(i,:) - [xsc,y(i),delta(i)];
        J(i,:) = J(i,:) - [xsc,y(i),delta(i)];
        Rot = [cos(alpha(i)), 0, sin(alpha(i));...
                  0,            1,  0          ;...
              -sin(alpha(i)), 0, cos(alpha(i))];
        A(i,:) = (Rot*A(i,:)')';
        B(i,:) = (Rot*B(i,:)')';
        C(i,:) = (Rot*C(i,:)')';
        D(i,:) = (Rot*D(i,:)')';
        E(i,:) = (Rot*E(i,:)')';
        F(i,:) = (Rot*F(i,:)')';
        G(i,:) = (Rot*G(i,:)')';
        H(i,:) = (Rot*H(i,:)')';
        I(i,:) = (Rot*I(i,:)')';
        J(i,:) = (Rot*J(i,:)')';
        A(i,:) = A(i,:) + [xsc,y(i),delta(i)];
        B(i,:) = B(i,:) + [xsc,y(i),delta(i)];
        C(i,:) = C(i,:) + [xsc,y(i),delta(i)];
        D(i,:) = D(i,:) + [xsc,y(i),delta(i)];
        E(i,:) = E(i,:) + [xsc,y(i),delta(i)];
        F(i,:) = F(i,:) + [xsc,y(i),delta(i)];
        G(i,:) = G(i,:) + [xsc,y(i),delta(i)];
        H(i,:) = H(i,:) + [xsc,y(i),delta(i)];
        I(i,:) = I(i,:) + [xsc,y(i),delta(i)];
        J(i,:) = J(i,:) + [xsc,y(i),delta(i)];
%         plot3( [A(i,1),B(i,1)],...
%                [A(i,2),B(i,2)],...
%                [A(i,3),B(i,3)],'b-');
%         plot3( [A(i,1),C(i,1)],...
%                [A(i,2),C(i,2)],...
%                [A(i,3),C(i,3)],'b-');
        plot3( [B(i,1),C(i,1)],...
               [B(i,2),C(i,2)],...
               [B(i,3),C(i,3)],'b-');
        plot3( [D(i,1),E(i,1)],...
               [D(i,2),E(i,2)],...
               [D(i,3),E(i,3)],'b-');
        plot3( [D(i,1),F(i,1)],...
               [D(i,2),F(i,2)],...
               [D(i,3),F(i,3)],'b-');
        plot3( [D(i,1),F(i,1)],...
               [D(i,2),F(i,2)],...
               [D(i,3),F(i,3)],'b-');
        plot3( [G(i,1),F(i,1)],...
               [G(i,2),F(i,2)],...
               [G(i,3),F(i,3)],'b-');
        plot3( [G(i,1),E(i,1)],...
               [G(i,2),E(i,2)],...
               [G(i,3),E(i,3)],'b-');
        plot3( [H(i,1),I(i,1)],...
               [H(i,2),I(i,2)],...
               [H(i,3),I(i,3)],'b-');
        plot3( [H(i,1),J(i,1)],...
               [H(i,2),J(i,2)],...
               [H(i,3),J(i,3)],'b-');
        plot3( [J(i,1),I(i,1)],...
               [J(i,2),I(i,2)],...
               [J(i,3),I(i,3)],'b-');
           
    end
    plot3(A(:,1),A(:,2),A(:,3),'b-')
    plot3(B(:,1),B(:,2),B(:,3),'b-')
    plot3(C(:,1),C(:,2),C(:,3),'b-')
    plot3(D(:,1),D(:,2),D(:,3),'b-')
    plot3(E(:,1),E(:,2),E(:,3),'b-')
    plot3(F(:,1),F(:,2),F(:,3),'b-')
    plot3(G(:,1),G(:,2),G(:,3),'b-')
    plot3(H(:,1),H(:,2),H(:,3),'b-')
    plot3(I(:,1),I(:,2),I(:,3),'b-')
    plot3(J(:,1),J(:,2),J(:,3),'b-')
    % For each rib
    for i=1:length(y_ribs)-1
        % get the parameters associated to this rib
        myAlpha = alpha(Tn(i,1));
        myY     = y(Tn(i,1));
        myDelta = delta(Tn(i,1));
        
        % for each point in the rib
        for j=1:N
           % Get its coordinates from the origin
           NACA_point_sup(j,:) = [x_for_NACAs(j),...
                                  myY,...
                                  myNACA(x_for_NACAs(j))+myDelta];
           NACA_point_inf(j,:) = [x_for_NACAs(j),...
                                  myY,...
                                  -myNACA(x_for_NACAs(j))+myDelta];
           % Do a translation to the xsc, as we will rotate around it
           NACA_point_sup(j,:) = NACA_point_sup(j,:) - [xsc,myY,myDelta];
           NACA_point_inf(j,:) = NACA_point_inf(j,:) - [xsc,myY,myDelta];                  
           % My rotation matrix
           Rot = [cos(myAlpha), 0, sin(myAlpha);...
                  0,            1,  0          ;...
                  -sin(myAlpha), 0, cos(myAlpha)];
           % Rotate the points
           NACA_point_sup(j,:) = (Rot*NACA_point_sup(j,:)')';
           NACA_point_inf(j,:) = (Rot*NACA_point_inf(j,:)')';
           % Translate it back to the global coordinates
           NACA_point_sup(j,:) = NACA_point_sup(j,:) + [xsc,myY,myDelta];
           NACA_point_inf(j,:) = NACA_point_inf(j,:) + [xsc,myY,myDelta];
           NACA_point_sup1(j,:) = NACA_point_sup(j,:);
           NACA_point_sup2(j,:) = NACA_point_sup(j,:);
           NACA_point_inf1(j,:) = NACA_point_inf(j,:);
           NACA_point_inf2(j,:) = NACA_point_inf(j,:);
           NACA_point_sup1(j,2) = NACA_point_sup(j,2)-0.5*thickness_rib;
           NACA_point_inf1(j,2) = NACA_point_inf(j,2)-0.5*thickness_rib;
           NACA_point_sup2(j,2) = NACA_point_sup(j,2)+0.5*thickness_rib;
           NACA_point_inf2(j,2) = NACA_point_inf(j,2)+0.5*thickness_rib;
        end
        plot3(NACA_point_sup1(:,1),...
              NACA_point_sup1(:,2),...
              NACA_point_sup1(:,3),'k-');
        plot3(NACA_point_sup2(:,1),...
              NACA_point_sup2(:,2),...
              NACA_point_sup2(:,3),'k-');
        plot3(NACA_point_inf1(:,1),...
              NACA_point_inf1(:,2),...
              NACA_point_inf1(:,3),'k-');
        plot3(NACA_point_inf2(:,1),...
              NACA_point_inf2(:,2),...
              NACA_point_inf2(:,3),'k-');
    end
    axis equal
   
end
