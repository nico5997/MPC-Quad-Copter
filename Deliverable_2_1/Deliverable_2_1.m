clear all
close all
clc

quad = Quad();

% Generate Trimmed and Linearized system
[xs,us] = quad.trim(); % Compute steady?state for which 0 = f(xs,us)
sys = quad.linearize(xs, us); % Linearize the nonlinear model

% Matrices before transformation
quad = Quad();
A_1 = sys.A;
B_1 = sys.B;
C_1 = sys.C;
D_1 = sys.D;

% Apply the system transformation
sys_transformed = sys * inv(quad.T); % New system is A * x + B * inv(T) * v

% Matrices after tranformation

A_2 = sys_transformed.A;
B_2 = sys_transformed.B;
C_2 = sys_transformed.C;
D_2 = sys_transformed.D;

% Dynamics change test

A_1 == A_2
B_1 == B_2

% From this test, we can see that effectively the matrix A did not change while the
% matrix B did change. 

% Controllability test

Mc_F = ctrb(A_2,B_2(:,1));
Mc_F(Mc_F < 1e-15 & Mc_F > -1e-15) = 0;

Mc_M_alpha = ctrb(A_2,B_2(:,2));
Mc_M_alpha(Mc_M_alpha < 1e-15 & Mc_M_alpha > -1e-15) = 0;

Mc_M_beta = ctrb(A_2,B_2(:,3));
Mc_M_beta(Mc_M_beta < 1e-15 & Mc_M_beta > -1e-15) = 0;

Mc_M_gamma = ctrb(A_2,B_2(:,4));
Mc_M_gamma(Mc_M_gamma < 1e-15 & Mc_M_gamma > -1e-15) = 0;

% By calculating the controllability matrix for each input, we can deduce
% the variables that can be control by each input. For instance, the matrix 
% Mc_F shows non-zeros entity in in line 9 and 12, which corresponds to the
% variables z and derivative of z. Thus, if we would ignore the other
% variables the matrix would be full rank and the system controllable. F
% can only control z and its derivative.
% The same criteria can be observ for M_alpha with alpha, alpha_p, y, y_p.
% The same criteria can be observ for M_beta with beta, beta_p, x, x_p.
% The same criteria can be observ for M_gamma with gamma, gamma_p.

% Compute the 4 different systems
[sys_x, sys_y, sys_z, sys_yaw] = quad.decompose(sys, xs, us);


