clear
close all

% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

% 1.1 - Input data
% Select mass to perform the loop for finding the shear center
mass_test = 1;
% location of colloction points for masses
mass_location = -0.1:0.02:0.2;
% if desired the plot for finding the x_sc
lets_plot = 0;
% slices to calculate the shear center
slices_sc = 8:14; % haaaardcoooode

%% STEP 1 - Finding the average shear center position x_sc

% 1.2 - Loop over the collocation points for masses and obtain the torsion
mass_points = length(mass_location);
theta_last = zeros(1,mass_points);
for i=1:mass_points
   % 1.2.1 - Set loads to zeros
   fload = zeros(1,mass_points);
   % 1.2.2 - Collocate the mass at the point desired
   fload(i) = mass_test;
   % 1.2.3 - Solve the system for the imposed BC
   displ = PazyWingLoad(fload,false);
   % 1.2.4 - Apply formula in slide 8 to find the torsion at the last
   %         slide
   theta_last(i) = (displ(end-1,1)-displ(end,1))/60;
end

% 1.3 - Detect the hange in sign of theta
sign = find(theta_last(2:end).*theta_last(1:end-1)<0);
x_sc = interp1([theta_last(sign),theta_last(sign+1)],...
    [mass_location(sign),mass_location(sign+1)],0);
if x_sc<0 || x_sc>90, error('A la parra x_sc = %f',x_sc); end

% 1.4 - Do a not much beautyful plot 
if lets_plot    
    figure
    hold on
    grid on
    plot(mass_location,theta_last,'LineWidth',1.5)
    xlabel('Mass location [mm]')
    ylabel('Elastic twist of last section \theta');
    xline(x_sc)
    text(x_sc+mean(mass_location)*0.2,max(theta_last)*0.75,'x_{sc}')
end

%% STEP 2 - Evaluation of twist and displacement for each y for each mass

% 2.1 - loop for the different possible mass positions
for i=1:mass_points
   % 2.1.1 - Set loads to zeros
   fload = zeros(1,mass_points);
   % 2.1.2 - Collocate the mass at the point desired
   fload(i) = mass_test;
   % 2.1.3 - Solve the system for the imposed BC
   displ = PazyWingLoad(fload,false);
   % 2.1.4 - Apply formula in slide 8 to find the torsion at the last
   %         slice
   theta(i,:) = (displ(1:2:end-1,1)-displ(2:2:end,1))/60;
end

% 2.2 - loop for each slice
for j=1:size(theta,2)
    % 2.2.1 - Calculate the sign change
    sign = find(theta(2:end,j).*theta(1:end-1,j)<0);
    % WARNING!!! Slices with multiple crosses to 0 are skipped and their
    % x_sc is set to zero!
    if length(sign)==1
        x_sc(j) = interp1([theta(sign,j),theta(sign+1,j)],...
        [mass_location(sign),mass_location(sign+1)],0);
    end
end

% 2.3 - Check the results
if min(x_sc)<0 || max(x_sc)>90
    error('A la parra x_sc <min = %f, max =5f>',min(x_sc),max(x_sc)); 
end

% 2.4 - Plot the x_sc at each slice for a first check
if 1
    figure
    hold on
   plot(x_sc,'LineWidth',1.5)
   grid on
   xlabel('i slice')
   ylabel('x_sc')
end

estimated.x_sc = mean(x_sc(slices_sc));

%% STEP 3 - Computation of the effective torsional stifness

% 3.1 - Make a vector with the span position of each slice
y_slices = displ(1:2:end,3);
if length(y_slices)~=size(theta,2)
    error('theta does not correspont to y_slices, what are you doing?')
end

for i=1:mass_points
    [a(i,:),S1(i)] = polyfit(y_slices,(theta(i,:))',1);  
    dThetady(i) = a(i,1);
    for j=1:length(y_slices)
        GJ(i,j) = (mass_test*9.81*(mass_location(i)-x_sc(j)))/dThetady(i);
%         if GJ(i,j)<0, error('OMG GJ(%d,%d)=%f<00',i,j,GJ(i,j)); end
    end
end
figure
surf(GJ);
estimated.GJ = mean(mean(GJ(:,[1:6,12:end])))/1e3; % hardcode 1e3

%% STEP 4 - Computation of the bending stifness
for i=1:mass_points
    fload = zeros(1,mass_points);
     fload(i) = mass_test;
    displ = PazyWingLoad(fload,false);
    count = 1;
    for j=1:2:size(displ,1)-1      
        w_sc(count) = interp1([displ(j,2),displ(j+1,2)],...
        [displ(j,1),displ(j+1,1)],x_sc(count));
        count = count +1;
    end
    [b(i,:),S1(i)] = polyfit(y_slices,w_sc,3);  
    d3wscdy3(i) = b(i,1);

    EI(i) = (mass_test*9.81)/(6*d3wscdy3(i));
%     if EI(i)<0, error('WTF EI negative in mass point i=%d',i); end
end
figure
plot(EI,'LineWidth',1.5)
estimated.EI =  mean(EI(slices_sc));  

estimated

save('stiffness/estimated.mat','estimated')
