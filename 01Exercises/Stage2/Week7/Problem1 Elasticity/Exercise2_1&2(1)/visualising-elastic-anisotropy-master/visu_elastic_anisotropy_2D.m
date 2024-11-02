clc;
clear all;


C = [     0.0905    0.0641    0.0001;
    0.0641    0.0904    0.0001;
    0.0001    0.0001    0.0095];
% C = [   0.5553    0.1776    0.0000
%         0.1776    0.5569   -0.0000
%         0.0000   -0.0000    0.0899];
C11=C(1,1);C22=C(2,2);C12=C(1,2);C66=C(3,3);C16=C(1,3);C26=C(2,3);
detC=C11*C22-C12^2;
theta = linspace(0,2*pi,360);
c=cos(theta);
s=sin(theta);

% Young's modulus
E =1./(C11*c.^4+C22*s.^4+(2*C12+C66).*c.^2.*s.^2+2*C16*c.^3.*s+2*C26.*s.^3.*c);
x_E=E.*c;
y_E=E.*s;
% % Bulk modulus
% K         = zeros(1,n);
% % Normal vector
% nu        = zeros(1,n);

% Shear modulus
G=4*(C11+C22-2*C12)*c.^2.*s.^2+C66*(c.^4+s.^4-2*s.^2.*c.^2);
G= 1./G;

%% Plot of Young's modulus and Bulk modulus
% Young's modulus    
a = figure(1);
% Plot of the surface in cartesian coordiantes. The fourth argument 
% "E" respectively "K" is for the correct color.
polarplot(theta, G);             
