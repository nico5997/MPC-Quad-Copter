clear all, close all, clc

%% Part 3

Ts = 1/5;
quad = Quad(Ts);
[xs, us] = quad.trim();
sys = quad.linearize(xs, us);
[sys_x, sys_y, sys_z, sys_yaw] = quad.decompose(sys, xs, us);

global tmp;

tmp = 0;

% Design MPC controller


N = 100;
Work_consumption_x = 0;
Work_consumption_y = 0;
Work_consumption_z = 0;
Work_consumption_yaw = 0;


x = zeros(4,N);
ux = zeros(1,N);
ux_tmp = zeros(1,N);
x(:,1) = [0 0 0 2]';
y = zeros(4,N);
uy = zeros(1,N);
uy_tmp = zeros(1,N);
y(:,1) = [0 0 0 2]';
z = zeros(2,N);
uz = zeros(1,N);
uz_tmp = zeros(1,N);
z(:,1) = [0 2]';
yaw = zeros(2,N);
uyaw = zeros(1,N);
uyaw_tmp = zeros(1,N);
yaw(:,1) = [0 pi/4]';

while (tmp < 2)
    mpc_x = MPC_Control_x(sys_x, Ts);
    mpc_y = MPC_Control_y(sys_y, Ts);
    mpc_z = MPC_Control_z(sys_z, Ts);
    mpc_yaw = MPC_Control_yaw(sys_yaw, Ts);
    for i = 1:N-1

        % Extract the optimal input
        ux(:,i) = mpc_x.get_u(x(:,i));
        uy(:,i) = mpc_y.get_u(y(:,i));
        uz(:,i) = mpc_z.get_u(z(:,i));
        uyaw(:,i) = mpc_yaw.get_u(yaw(:,i));

        if (tmp == 0)
           ux_tmp(:,i) = mpc_x.get_u(x(:,i));
           uy_tmp(:,i) = mpc_y.get_u(y(:,i));
           uz_tmp(:,i) = mpc_z.get_u(z(:,i));
           uyaw_tmp(:,i) = mpc_yaw.get_u(yaw(:,i)); 
        end

        if (i ~= 1)
            Work_consumption_x = Work_consumption_x + abs(ux(1,i)*(x(4,i)-x(4,i-1)));
            Work_consumption_y = Work_consumption_y + abs(uy(1,i)*(y(4,i)-y(4,i-1)));
            Work_consumption_z = Work_consumption_z + abs(uz(1,i)*(z(2,i)-z(2,i-1)));
            Work_consumption_yaw = Work_consumption_yaw + abs(uyaw(1,i)*(yaw(2,i)-yaw(2,i-1)));
        end
    i
        % Apply the optimal input to the system
        x(:,i+1) = mpc_x.A*x(:,i) + mpc_x.B*ux(:,i);
        y(:,i+1) = mpc_y.A*y(:,i) + mpc_y.B*uy(:,i);
        z(:,i+1) = mpc_z.A*z(:,i) + mpc_z.B*uz(:,i);
        yaw(:,i+1) = mpc_yaw.A*yaw(:,i) + mpc_yaw.B*uyaw(:,i);
    end
    tmp = tmp + 1;
end

% Plots of each dimentions

figure
grid on

subplot(411);
plot(Ts*(1:N),x(4,:));
ylabel('Position x (m)');
xlabel('time [s]');
title('Position x of Quadcopter with respect to time')

subplot(412);
plot(Ts*(1:N),y(4,:));
ylabel('Position y [m]');
xlabel('time [s]');
title('Position y of Quadcopter with respect to time')

subplot(413);
plot(Ts*(1:N),z(2,:));
ylabel('Position z [m]');
xlabel('time [s]');
title('Position z of Quadcopter with respect t0 time')

subplot(414);
plot(Ts*(1:N),yaw(2,:));
ylabel('Angle yaw [rad]');
xlabel('time [s]');
title('Angle yaw of Quadcopter with respect to time')

S_x = stepinfo(x(4,:),Ts*(1:N),0,'SettlingTimeThreshold',0.02);
disp('Settling time of x is:')
S_x.SettlingTime
S_y = stepinfo(y(4,:),Ts*(1:N),0,'SettlingTimeThreshold',0.02);
disp('Settling time of y is:')
S_y.SettlingTime
S_z = stepinfo(z(2,:),Ts*(1:N),0,'SettlingTimeThreshold',0.02);
disp('Settling time of z is:')
S_z.SettlingTime
S_yaw = stepinfo(yaw(2,:),Ts*(1:N),0,'SettlingTimeThreshold',0.02);
disp('Settling time of yaw is:')
S_yaw.SettlingTime

% Plots of contraints

% Constraints over the input
figure
grid on

subplot(411);
plot(Ts*(1:N),0.3*ones(1,length(ux(1,:))),'r--');
hold on
plot(Ts*(1:N),-0.3*ones(1,length(ux(1,:))),'r--');
plot(Ts*(1:N),ux(1,:));
ylabel('Moment M\_beta [Nm]');
xlabel('time [s]');
axis([1/5 10 -0.4 0.4]);
title('Input moment M\_beta of Quadcopter with respect to time')
legend('Contraints')

subplot(412);
plot(Ts*(1:N),0.3*ones(1,length(ux(1,:))),'r--');
hold on
plot(Ts*(1:N),-0.3*ones(1,length(ux(1,:))),'r--');
plot(Ts*(1:N),uy(1,:));
ylabel('Moment M\_alpha [Nm]');
xlabel('time [s]');
axis([1/5 10 -0.4 0.4]);
title('Input moment M\_alpha of Quadcopter with respect to time')
legend('Contraints')

subplot(413);
plot(Ts*(1:N),0.3*ones(1,length(ux(1,:))),'r--');
hold on
plot(Ts*(1:N),-0.2*ones(1,length(ux(1,:))),'r--');
plot(Ts*(1:N),uz(1,:));
ylabel('Force F [N]');
xlabel('time [s]');
axis([1/5 10 -0.4 0.4]);
title('Input force F of Quadcopter with respect to time')
legend('Contraints')

subplot(414);
plot(Ts*(1:N),0.2*ones(1,length(ux(1,:))),'r--');
hold on
plot(Ts*(1:N),-0.2*ones(1,length(ux(1,:))),'r--');
plot(Ts*(1:N),uyaw(1,:));
ylabel('Moment M\_gamma [Nm]');
xlabel('time [s]');
axis([1/5 10 -0.4 0.4]);
title('Input moment M\_gamma of Quadcopter with respect to time')
legend('Contraints')

% Comparison between u and u_optimized

figure
grid on

subplot(411);
plot(Ts*(1:N),ux_tmp(1,:));
hold on
plot(Ts*(1:N),ux(1,:));
ylabel('Moment M\_beta [Nm]');
xlabel('time [s]');
axis([1/5 10 -0.4 0.4]);
title('Input moment M\_beta of Quadcopter with respect to time')
legend('input','optimized input')

subplot(412);
plot(Ts*(1:N),uy_tmp(1,:));
hold on
plot(Ts*(1:N),uy(1,:));
ylabel('Moment M\_alpha [Nm]');
xlabel('time [s]');
axis([1/5 10 -0.4 0.4]);
title('Input moment M\_alpha of Quadcopter with respect to time')
legend('input','optimized input')

subplot(413);
plot(Ts*(1:N),uz_tmp(1,:));
hold on
plot(Ts*(1:N),uz(1,:));
ylabel('Force F [N]');
xlabel('time [s]');
axis([1/5 10 -0.4 0.4]);
title('Input force F of Quadcopter with respect to time')
legend('input','optimized input')

subplot(414);
plot(Ts*(1:N),uyaw_tmp(1,:));
hold on
plot(Ts*(1:N),uyaw(1,:));
ylabel('Moment M\_gamma [Nm]');
xlabel('time [s]');
axis([1/5 10 -0.4 0.4]);
title('Input moment M\_gamma of Quadcopter with respect to time')
legend('input','optimized input')

% constraints over the states

figure
grid on

subplot(211);
plot(Ts*(1:N),0.035*ones(1,length(ux(1,:))),'r--');
hold on
plot(Ts*(1:N),-0.035*ones(1,length(ux(1,:))),'r--');
plot(Ts*(1:N),x(2,:));
ylabel('Angle beta [rad]');
xlabel('time [s]');
axis([1/5 10 -0.05 0.05]);
title('Angles beta of Quadcopter with respect to time')
legend('Contraints')

subplot(212);
plot(Ts*(1:N),0.035*ones(1,length(ux(1,:))),'r--');
hold on
plot(Ts*(1:N),-0.035*ones(1,length(ux(1,:))),'r--');
plot(Ts*(1:N),y(2,:));
ylabel('Angle alpha [rad]');
xlabel('time [s]');
axis([1/5 10 -0.05 0.05]);
title('Angles alpha of Quadcopter with respect to time')
legend('Contraints')