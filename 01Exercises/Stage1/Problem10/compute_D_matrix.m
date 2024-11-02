
% Define the material properties
E = 100;  % Young's modulus in Pascals (e.g., 210 GPa for steel)
mu = 0.3;   % Poisson's ratio (e.g., 0.3 for typical metals)

% Construct the elasticity matrix D for plane stress condition
D = (E / (1 - mu^2)) * [1, mu, 0;
                        mu, 1, 0;
                        0, 0, (1 - mu) / 2];

% Display the result
disp('The D matrix is:');
disp(D);
