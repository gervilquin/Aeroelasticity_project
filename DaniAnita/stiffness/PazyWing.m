%% Virtual test Pazy wing

clear
close all

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

% Test is performed and vertical displacements are recorded in displ(:,1).
% displ(:,2) are the x-coordinates of each sensor (chordwise direction).
% displ(:,3) are the y-coordinates of each sensor (spanwise direction).
displ = PazyWingLoad(fload,true);

% Note: to prevent plot from showing set second input in the 'PazyWingLoad'
% funciton to 'false'.