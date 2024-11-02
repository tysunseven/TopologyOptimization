% Define constants
a = 0.5;
b = 0.5;

% Define symbolic variables
syms x y real

% Define the shape functions for a four-node rectangular element
N1 = (1/4) * (1 + x/a) * (1 + y/b);
N2 = (1/4) * (1 - x/a) * (1 + y/b);
N3 = (1/4) * (1 - x/a) * (1 - y/b);
N4 = (1/4) * (1 + x/a) * (1 - y/b);

% Construct the shape function matrix N using the defined shape functions
N = [N1 0 N2 0 N3 0 N4 0; 
     0 N1 0 N2 0 N3 0 N4];

% Calculate partial derivatives of the shape functions
dN1_dx = diff(N1, x); dN1_dy = diff(N1, y);
dN2_dx = diff(N2, x); dN2_dy = diff(N2, y);
dN3_dx = diff(N3, x); dN3_dy = diff(N3, y);
dN4_dx = diff(N4, x); dN4_dy = diff(N4, y);

% Substitute x=0 and y=0 into the partial derivatives
dN1_dx_val = subs(dN1_dx, {x, y}, {0, 0});
dN1_dy_val = subs(dN1_dy, {x, y}, {0, 0});
dN2_dx_val = subs(dN2_dx, {x, y}, {0, 0});
dN2_dy_val = subs(dN2_dy, {x, y}, {0, 0});
dN3_dx_val = subs(dN3_dx, {x, y}, {0, 0});
dN3_dy_val = subs(dN3_dy, {x, y}, {0, 0});
dN4_dx_val = subs(dN4_dx, {x, y}, {0, 0});
dN4_dy_val = subs(dN4_dy, {x, y}, {0, 0});

% Construct the B matrix using the evaluated partial derivatives
B = [dN1_dx_val 0 dN2_dx_val 0 dN3_dx_val 0 dN4_dx_val 0;
     0 dN1_dy_val 0 dN2_dy_val 0 dN3_dy_val 0 dN4_dy_val;
     dN1_dy_val dN1_dx_val dN2_dy_val dN2_dx_val dN3_dy_val dN3_dx_val dN4_dy_val dN4_dx_val];

% Define the material properties
E = 100;  % Young's modulus in Pascals (e.g., 210 GPa for steel)
mu = 0.3;   % Poisson's ratio (e.g., 0.3 for typical metals)

% Construct the elasticity matrix D for plane stress condition
D = (E / (1 - mu^2)) * [1, mu, 0;
                        mu, 1, 0;
                        0, 0, (1 - mu) / 2];

S=D*B;

% Display the final B matrix
disp('The S matrix is:');
disp(S);
