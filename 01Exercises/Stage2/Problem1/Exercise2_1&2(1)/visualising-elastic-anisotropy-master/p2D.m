clc;
clear all;

%%111.jpg
C = 1e-3*[ 0.4333 -0.1668 0.0980;
           -0.1668 0.4609 -0.1004;
            0.0980 -0.1004 0.0590 ];
%%222.jpg
% C= 1.0e-03 *[ 0.6358   -0.3889   -0.0017;
%              -0.3889    0.6263    0.0016;
%              -0.0017    0.0016    0.0224];

%micro_1.png
%  C=[  0.2866    0.1056   -0.0000
%       0.1056    0.2858   -0.0000
%      -0.0000   -0.0000    0.0127  ];

%micro_2.png
% C=[ 0.5553    0.1776    0.0000
%     0.1776    0.5569   -0.0000
%     0.0000   -0.0000    0.0899];

%micro_3.png
% C=[ 0.1371    0.1114   -0.0000
%     0.1114    0.1354   -0.0000
%    -0.0000   -0.0000    0.0216  ];
C11=C(1,1);C22=C(2,2);C12=C(1,2);C66=C(3,3);C16=C(1,3);C26=C(2,3);
detC=C11*C22-C12^2;
theta = linspace(0,2*pi,360);
c=cos(theta);
s=sin(theta);

% Young's modulus
E =1./((C22*c.^4-2*C12.*c.^2.*s.^2+C11*s.^4)/detC+c.^2.*s.^2/C66);
v =(C12*c.^4-(C11+C22-detC/C66)*c.^2.*s.^2+C12*s.^4)./(C22*c.^4-(2*C12-detC/C66)*c.^2.*s.^2+C11*s.^4);
vp_index=(v>0);
vn_index=(v<=0);

% Shear modulus
G=4*(C11+C22-2*C12)*c.^2.*s.^2+C66*(c.^4+s.^4-2*s.^2.*c.^2);
%G= 1./G;
G = Shear_2DM(C, theta);

%% Plot of Young's modulus and Bulk modulus
% Young's modulus    
a = figure(1);
% Plot of the surface in cartesian coordiantes. The fourth argument 
% "E" respectively "K" is for the correct color.
polarplot(theta, G);   
% polarplot(theta, E);