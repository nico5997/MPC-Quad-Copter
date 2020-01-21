classdef MPC_Control_z < MPC_Control
  properties
    A_bar, B_bar, C_bar % Augmented system for disturbance rejection    
    L                   % Estimator gain for disturbance rejection
  end
  
  methods
    function mpc = MPC_Control_z(sys, Ts)
      mpc = mpc@MPC_Control(sys, Ts);
      
      [mpc.A_bar, mpc.B_bar, mpc.C_bar, mpc.L] = mpc.setup_estimator();
    end
    
    % Design a YALMIP optimizer object that takes a steady-state state
    % and input (xs, us) and returns a control input
    function ctrl_opt = setup_controller(mpc)

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % INPUTS
      %   x(:,1) - initial state (estimate)
      %   d_est  - disturbance estimate
      %   xs, us - steady-state target
      % OUTPUTS
      %   u(:,1) - input to apply to the system
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      [n,m] = size(mpc.B);
      
      % Steady-state targets (Ignore this before Todo 3.3)
      xs = sdpvar(n, 1);
      us = sdpvar(m, 1);
      
      % Disturbance estimate (Ignore this before Part 5)
      d_est = sdpvar(1);

      % SET THE HORIZON HERE
      N = 10;
      
      % Predicted state and input trajectories
      x = sdpvar(n, N);
      u = sdpvar(m, N-1);
      

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE 

      % NOTE: The matrices mpc.A, mpc.B, mpc.C and mpc.D are 
      %       the DISCRETE-TIME MODEL of your system

      % SET THE PROBLEM CONSTRAINTS con AND THE OBJECTIVE obj HERE
      con = [];
      obj = 0;
      % SET PARAMETERS HERE
      Q = [20 0;
           0 800];
      R = 1*eye(m);
      
      % CONSTRAINTS
      %  u in U = { u | Mu <= m }
      M = [1;-1]; m = [0.3; 0.2];
      % x in X = { x | Fx <= f }
      F = [0 0; 0 0]; f = [0; 0];
      
      % Compute LQR controller for unconstrained system
      [K,Qf,~] = dlqr(mpc.A,mpc.B,Q,R);
      % MATLAB defines K as -K, so invert its signal
      K = -K; 
      eig(mpc.A + mpc.B*K)

      % Compute maximal invariant set
      Xf = polytope([F;M*K],[f;m]);
      Acl = [mpc.A + mpc.B*K];
      while 1
          prevXf = Xf;
          [T,t] = double(Xf);
          preXf = polytope(T*Acl,t);
          Xf = intersect(Xf, preXf);
          if isequal(prevXf, Xf)
              break
          end
      end
      [Ff,ff] = double(Xf);
      
      con = (x(:,2) == mpc.A*x(:,1) + mpc.B*u(:,1) + mpc.B*d_est) + (M*u(:,1) <= m);
      obj = (u(:,1)-us)'*R*(u(:,1)-us);
      for i = 2:N-1
          con = con + (x(:,i+1) == mpc.A*x(:,i) + mpc.B*u(:,i) + mpc.B*d_est);
          con = con + (M*u(:,i) <= m);
          obj = obj + (x(:,i)-xs)'*Q*(x(:,i)-xs) + (u(:,i)-us)'*R*(u(:,i)-us);
      end
      %con = con + (Ff*x(:,N) <= ff);
      obj = obj + (x(:,N)-xs)'*Qf*(x(:,N)-xs);

      
      % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE 
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      
      ctrl_opt = optimizer(con, obj, sdpsettings('solver','gurobi'), ...
        {x(:,1), xs, us, d_est}, u(:,1));
    end
    
    
    % Design a YALMIP optimizer object that takes a position reference
    % and returns a feasible steady-state state and input (xs, us)
    function target_opt = setup_steady_state_target(mpc)
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % INPUTS
      %   ref    - reference to track
      %   d_est  - disturbance estimate
      % OUTPUTS
      %   xs, us - steady-state target
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      % Steady-state targets
      n = size(mpc.A,1);
      xs = sdpvar(n, 1);
      us = sdpvar;
      
      % Reference position (Ignore this before Todo 3.3)
      ref = sdpvar;
      
      % Disturbance estimate (Ignore this before Part 5)
      d_est = sdpvar(1);
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE 
      % You can use the matrices mpc.A, mpc.B, mpc.C and mpc.D
      con = [];
      obj = 0;
      % SET PARAMETERS HERE
      Qs = 1;
      Rs = 10;
            
      % CONSTRAINTS
      %  u in U = { u | Mu <= m }
      M = [1;-1]; m = [0.3; 0.2];
      % x in X = { x | Fx <= f }
      F = [0 0; 0 0]; f = [0; 0];
      
      %con = (xs == mpc.A*xs + mpc.B*us);
      %con = con + (M*us <= m);
      %obj = (mpc.C*xs - ref)'*Qs*(mpc.C*xs-ref);
      
      con = ([eye(n)-mpc.A, -mpc.B; mpc.C, zeros(size(mpc.D))]*[xs;us] == [-mpc.B*d_est; ref - d_est]);
      con = con + (M*us <= m);
      obj = us'*Rs*us;
      
      % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE 
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      
      % Compute the steady-state target
      target_opt = optimizer(con, obj, sdpsettings('solver', 'gurobi'), {ref, d_est}, {xs, us});
    end
    
    
    % Compute augmented system and estimator gain for input disturbance rejection
    function [A_bar, B_bar, C_bar, L] = setup_estimator(mpc)
      
      %%% Design the matrices A_bar, B_bar, L, and C_bar
      %%% so that the estimate x_bar_next [ x_hat; disturbance_hat ]
      %%% converges to the correct state and constant input disturbance
      %%%   x_bar_next = A_bar * x_bar + B_bar * u + L * (C_bar * x_bar - y);
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE 
      % You can use the matrices mpc.A, mpc.B, mpc.C and mpc.D
      

      [n,m] = size(mpc.B);
      ny = size(mpc.C,2);
      
      A_bar = [mpc.A     ,  mpc.B    ;
               zeros(m,n),  eye(m,m)];
      B_bar = [mpc.B; zeros(1,m)];
      C_bar = [mpc.C, 1];
      
      
      O = obsv(A_bar,C_bar);
      O
      rank_obsv = rank(O) % full rank -> observable system
      
      L = -place(A_bar',C_bar',[0.2,0.3,0.1])';
      L
  
      
      % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE 
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

    
  end
end
