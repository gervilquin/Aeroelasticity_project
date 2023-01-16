%% Create structura mesh
%   input
% ns           -> number of nodes in each subset (int vector)
%   output
% y            -> global position of the nodes
% y_ribs       -> global polition of the ribs
% l            -> element length
% Authors: Ana María Díaz, Daniel Regener
% Advanced Aeroelasticity, ESEIAAT, 28 dec 2022


function [y, y_ribs, l] = compute_nodes(ns)
    
    % position of each rib
    y_ribs = 1e-3*[0. 38.4 76.8 115.2 153.6 192. 230.4 268.8 307.2 345.6 384. 422.4 460.8 499.2 537.6 541.6 545.6 550.];
    
    delta_ribs = zeros(1,length(y_ribs)-1);

    % distance between ribs
    for i=1:length(y_ribs)-1
        delta_ribs(i) = y_ribs(i+1) - y_ribs(i);
    end

    % position of each node
    y(1) = 0.;
    for e=1:size(ns,2)
        for i=2:ns(e)
            deltay = delta_ribs(e)/(ns(e)-1);
            yloc = (i-1)*deltay;
            y = [y, y_ribs(e)+yloc];
        end
    end
   
   % length of each element
   l = zeros(1,length(y)-1); 
   for i = 1:length(y)-1
       l(i)=y(i+1)-y(i);
   end
    
end