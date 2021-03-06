clear all, close all, clc

%% Part 5

BIAS = -0.1;
Ts = 1/5;
quad = Quad_modif(Ts); % to plot the modified path of z
%quad = Quad(Ts); % to plot the MPC path of z
[xs, us] = quad.trim();
sys = quad.linearize(xs, us);
[sys_x, sys_y, sys_z, sys_yaw] = quad.decompose(sys, xs, us);
mpc_x = MPC_Control_x(sys_x, Ts);
mpc_y = MPC_Control_y(sys_y, Ts);
mpc_z = MPC_Control_z(sys_z, Ts);
mpc_yaw = MPC_Control_yaw(sys_yaw, Ts);

sim = quad.sim(mpc_x, mpc_y, mpc_z, mpc_yaw, BIAS);
quad.plot(sim);
