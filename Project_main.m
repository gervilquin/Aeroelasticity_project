%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aeroelastic model of a Pazy wing
% Authors: Gerard Villalta & Antoni Alberto
% Date: 18/11/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc, clear, close all
addpath('functions\')

%% 1. Initialisation 

% Elmements definition
%                   Rib     Node j      yj  
y_sections =       [0,       1,           0.0
                    1,       2,           38.4
                    2,       3,           76.8
                    3,       4,           115.2
                    4,       5,           153.6
                    5,       6,           192.0
                    6,       7,           230.4
                    7,       8,           268.8
                    8,       9,           307.2
                    9,       10,          345.6
                    10,      11,          384.0
                    11,      12,          422.4
                    12,      13,          460.8
                    13,      14,          499.2
                    14,      15,          537.6
                    15,      16,          541.6
                    16,      17,          545.6
                    0,       18,          550.0 ];

y_sections = y_sections(:,2:3);

nsec = length(y_sections(:,2))-1; % number of panels
neset = 10;%<<<<<<<<<< Input 
nel = nsec*neset; % total number of elements

% Centers
x_sc = 0.43;
x_ac = 1/4;
x_col = 3/4;

% Geometry
chord = 100; %mm

% Material properties
EI = 6500;
GJ = 5500;

%% 2. Structural modelling

% Define the nodal coordinates
y_nodal = compute_x_nodal(y_sections(:,2),neset);
n = length(y_nodal);

% define the connectivity matrixs
Tn = ConnectivityElements(n-1);
Ts = ConnectivitySubsets(Tn,nsec-1,neset);

% Define stifness matrix
K = zeros(3*n);
for e=1:neset*nsec
    l = y_nodal(Tn)
end

%% 3. Aerodynamics modelling

% Define the chord vector
chord_y = ones(1,length(y_nodal))*chord;

% Compute area of the elements
S_ = compute_element_surface(chord_y,y_nodal);

% Compute chord/4 point and the colocation point span wise location
[ac_x,ac_y] = compute_aero_point(y_nodal,chord_y,x_ac);
[col_x,col_y] = compute_aero_point(y_nodal,chord_y,x_col);


%% 4. Aeroelastic linear coupling

%% 5. Aeroelastic solver

