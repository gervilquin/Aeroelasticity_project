%% Virtual test Pazy wing
clc;clear;close all

%% 1. Input data and initialization
y_section = [19.2; 57.6; 96.0; 134.4; 172.8; 211.2; 249.6; ...
            288.0; 326.4; 364.8; 403.2; 441.6; 480.0; 518.4];
n_section = length(y_section);
mass_pos = -0.1:0.02:0.2;
mass_slot = 1:1:16;

N_shear = 10;  % Experiment repetitions for each section (shear)
tip_torsion = zeros(1,16);


g = 9.81; 
N_rigid= 10; % Experiment repetitions for each section (rigidity)


theta = zeros(length(mass_slot),length(y_section));
wsc = zeros(length(mass_slot),length(y_section));
E_n = zeros(length(mass_slot),N_rigid);
G_n = zeros(length(mass_slot),N_rigid);

E = zeros(length(mass_slot),1);
G = zeros(length(mass_slot),1);

%% 2. Shear center analysis
Xsc_y = zeros(1,n_section);
for i = 1:n_section % Loop to compute effects on all sections
X_sc = zeros(1,N_shear);
for j=1:N_shear           % Loop for each experiment repetition
notfound = true;
k = 0;
while(notfound)
k = k+1;
% Choose masses (in kg) to add at each slot in the wing tip rod
fload = [
    0      % x01 = -0.10 m
    0      % x02 = -0.08 m
    0      % x03 = -0.06 m
    0      % x04 = -0.04 m
    0      % x05 = -0.02 m
    0      % x06 =  0.00 m
    0      % x07 =  0.02 m
    0      % x08 =  0.04 m
    0      % x09 =  0.06 m
    0      % x10 =  0.08 m
    0      % x11 =  0.10 m
    0      % x12 =  0.12 m
    0      % x13 =  0.14 m
    0      % x14 =  0.16 m
    0      % x15 =  0.18 m
    0      % x16 =  0.20 m
];

% Load vector generator (1kg placed at the different positions)
% Generates a binary vector changing the position of the 1.
Combination = 2^k;
fload = dec2bin(Combination,16) == '1'; 
fload = double(flip(fload'));

% Test is performed and vertical displacements are recorded in displ(:,1).
% displ(:,2) are the x-coordinates of each sensor (chordwise direction).
% displ(:,3) are the y-coordinates of each sensor (spanwise direction).
displ = PazyWingLoad(fload,false);

% Note: to prevent plot from showing set second input in the 'PazyWingLoad'
% function to 'false'.

notcvg = true;
A = 0;
A_old = -1;
increment = 1;
tip_torsion(k) = (displ((2*i-1),1) - displ(2*i,1))/60;
    if k>=2 
        if tip_torsion(k)*tip_torsion(k-1)<0 %When change in sign is found, start the iterative balancing process
            notfound = false;
            tip_torsion_old = tip_torsion(k);
            while(notcvg)
                fload = zeros(16,1);
                fload(k) = 1;
                fload(k-1) = A;
                displ = PazyWingLoad(fload,false);
                tip_torsion_zero = (displ((2*i-1),1) - displ(2*i,1))/60;
                if abs(tip_torsion_zero)<10e-8 && abs(A_old-A)<10e-6
                    notcvg = false;
                else
                    if tip_torsion_zero*tip_torsion_old<0 % If torsion changes polarity, too much mass was added.
                        A = A - increment;
                        increment = increment/10;
                    else
                        A = A + increment;
                    end
                end
                A_old = A;
            end
        end
    end
end
X_sc(j) = mass_pos(k-1) + 0.02/(A+1);
end
Xsc_y(i) = mean(X_sc);
end
Xsc = Xsc_y(end);
%% 3. Flexural and torsional rigidity
for i=1:length(mass_slot) % Loop over mass slots
    for j = 1:N_rigid % Loop over experiment repetition
        F_load = zeros(1,length(mass_slot));
        F_load(i) = i; 
        displ = PazyWingLoad(F_load,false);
        for k=1:length(y_section)  
            theta(i,k) = (displ((2*k-1),1) - displ(2*k,1))/0.06;
            wsc(i,k) = (displ(2*k,1) + theta(i,k)*(displ(2*k,2)-Xsc));
        end
        fit_theta = polyfit(y_section/1000,theta(i,:),1);
        fit_wsc = polyfit(y_section/1000,wsc(i,:),3);
        E_n(i,j) = -unique(-F_load(i))*g / (6 * fit_wsc(1));
        G_n(i,j) = unique(-F_load(i))*g * (Xsc - mass_pos(i))/fit_theta(1);
    end
    E(i) = mean(E_n(i,:));
    G(i) = mean(G_n(i,:));
end
E_mean = mean(E);
G_mean = mean([G(1),G(end)]);

%% 4. Results presentation
% Shear center position
figure(1)
hold on
plot(y_section,Xsc_y)
ylabel('Shear center position (m)')
xlabel('Wing position (m)')
xlim([0 520])
yline(min(Xsc_y),'--','color','black')
text(5,min(Xsc_y)+0.0015,'Shear center = 0.043')
hold off

% Flexural rigidity
figure()
plot(E,mass_pos)
xlabel('Flexrual rigidity (EI) [Pa·m4]')
ylabel('Position of the mass [mm]')
xline(E_mean,'--','color','black')
ylim([-0.1 0.2])
t = text(E_mean + 0.001,0.01,'EI = 5.04');
set(t,'Rotation',90)

% Torsional rigidity
figure()
plot(G,mass_pos)
xlabel('Torsional rigidity (GJ) [N·m2/rad]')
ylabel('Position of the mass [mm]')
xline(G_mean,'--','color','black')
ylim([-0.1 0.2])
t = text(G_mean + 0.1,0.1,'GJ = 6.59');
set(t,'Rotation',90)

% Torsion displacements plot
figure()
hold on
for i=1:length(theta(:,1))
    plot(y_section,theta(i,:))
end
xlabel('y')
ylabel('\theta')
hold off

% Bending displacements plot
figure()
hold on
for i=1:length(wsc(:,1))
    plot(y_section,wsc(i,:))
end
xlabel('y')
ylabel('w_{sc}')

%% 5. Data saving
StructureVEx.EI = E_mean;
StructureVEx.GJ = G_mean;
StructureVEx.Xsc = Xsc;

save('StructureVEx.mat','StructureVEx')