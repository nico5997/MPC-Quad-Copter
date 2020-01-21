function [ctrl, traj] = ctrl_NMPC(quad)

import casadi.*

opti = casadi.Opti(); % Optimization problem

N = 50; % MPC horizon [SET THIS VARIABLE]

% decision variables
X = opti.variable(12,N+1); % state trajectory variables
U = opti.variable(4, N);   % control trajectory
Xs = opti.variable(12,1);  % steady state reference (full state representation)

X0 = opti.parameter(12,1); % initial state
REF = opti.parameter(4,1); % reference position [x,y,z,yaw]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% YOUR CODE HERE %%%%

% ----- Build Xs and Us -----
[~, Us] = quad.trim();              % Grab sterady inputs for stationay flight
Xs([1:5,7:9]) =  zeros(8,1);         % These states are zero when hovering
Xs([10:12,6]) = REF;                % Grab reference position [x,y,z,yaw]

% ----- Constraints -----
opti.subject_to(0 <= U <= 1.5);     % Constraints on the input (from Todo 1.2)

% ----- Initial conditions -----
opti.subject_to(X(:,1) == X0);      % Initial condition on X

% ----- Objectives -----
Q = diag([0,0,0,0,0,1,0,0,0,100,100,10000]);   % Weight for cost on states
R = 3*diag([1,1,1,1]);                         % Weight for cost on control signals

obj = (U(:,1)-Us)'*R*(U(:,1)-Us); 
for i = 2:N-1
    obj = obj + (X(:,i)-Xs)'*Q*(X(:,i)-Xs);
    obj = obj + (U(:,i)-Us)'*R*(U(:,i)-Us);
end
obj = obj + (X(:,N)-Xs)'*Q*(X(:,N)-Xs);
opti.minimize(obj)

% ----- System dynamics -----
h = 1/5;                            % Sampling period (same as in part 3)

for k = 1:N % loop over control intervals
    k1 = quad.f(X(:,k),        U(:,k));
    k2 = quad.f(X(:,k)+h/2*k1, U(:,k));
    k3 = quad.f(X(:,k)+h/2*k2, U(:,k));
    k4 = quad.f(X(:,k)+h*k3,   U(:,k));
    x_next = X(:,k) + h/6*(k1+2*k2+2*k3+k4);
    opti.subject_to(X(:,k+1) == x_next);        % RK4 discretization
    
%   ----- Alternatively use Euler discretization -----
    %x_next_euler = X(:,k)+h*quad.f(X(:,k), U(:,k));
    %opti.subject_to(X(:,k+1) == x_next_euler)   % Euler discretization
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ctrl = @(x,ref) eval_ctrl(x, ref, opti, X0, REF, X, U);
end


function u = eval_ctrl(x, ref, opti, X0, REF, X, U)
% ---- Set the initial state and reference ----
opti.set_value(X0, x);
opti.set_value(REF, ref);

% ---- Setup solver NLP ----
ops = struct('ipopt', struct('print_level',0, 'tol', 1e-3), 'print_time', false);
opti.solver('ipopt', ops);

% ---- Solve the optimization problem ----
sol = opti.solve();
assert(sol.stats.success == 1, 'Error computing optimal input');

u = opti.value(U(:,1));

% Use the current solution to speed up the next optimization
opti.set_initial(sol.value_variables());
opti.set_initial(opti.lam_g, sol.value(opti.lam_g));
end