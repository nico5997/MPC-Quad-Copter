classdef MPC_Control_yaw < MPC_Control
  
  methods
    % Design a YALMIP optimizer object that takes a steady-state state
    % and input (xs, us) and returns a control input
    function ctrl_opt = setup_controller(mpc)

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % INPUTS
      %   x(:,1) - initial state (estimate)
      %   xs, us - steady-state target
      % OUTPUTS
      %   u(:,1) - input to apply to the system
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      [n,m] = size(mpc.B);
      
      % Steady-state targets (Ignore this before Todo 3.2)
      xs = sdpvar(n, 1);
      us = sdpvar(m, 1);
      
      % SET THE HORIZON HERE
      N = 10;
      
      % Predicted state and input trajectories
      x = sdpvar(n, N);
      u = sdpvar(m, N-1);
      

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE 

      % NOTE: The matrices mpc.A, mpc.B, mpc.C and mpc.D are 
      %       the DISCRETE-TIME MODEL of your system

      % WRITE THE CONSTRAINTS AND OBJECTIVE HERE
      con = [];
      obj = 0;
      % SET PARAMETERS HERE
      
      global tmp;
      
      if (tmp == 0)
         Q = [0 0;
              0 1];
         R = 1*eye(m); %20
      end
      
      if (tmp == 1)
          Q = [0 0;
              0 1];
          R = 20*eye(m); %20
      end
      
      
      
      % CONSTRAINTS
      %  u in U = { u | Mu <= m }
      M = [1;-1]; m = [0.2; 0.2];
      % x in X = { x | Fx <= f }
      F = [0 0; 0 0]; f = [0; 0];
      
      % Compute LQR controller for unconstrained system
      [K,Qf,~] = dlqr(mpc.A,mpc.B,Q,R);
      % MATLAB defines K as -K, so invert its signal
      K = -K; 

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
      
      if (tmp == 1)
      figure
      
      Xf.projection([1 2]).plot();
      xlabel('Angular speed gamma [rad/s]');
      ylabel('Angle gamma [rad]');
      title('Projection of Terminal invariant set: gamma\_dot vs. gamma')
      end
      
      con = (x(:,2) == mpc.A*(x(:,1)) + mpc.B*(u(:,1))) + (M*u(:,1) <= m);
      obj = (u(:,1))'*R*(u(:,1));
      for i = 2:N-1
          con = con + (x(:,i+1) == mpc.A*(x(:,i)) + mpc.B*(u(:,i)));
          con = con + (M*u(:,i) <= m);
          obj = obj + (x(:,i))'*Q*(x(:,i)) + (u(:,i))'*R*(u(:,i));
      end
      %con = con + (Ff*x(:,N) <= ff);
      obj = obj + (x(:,N))'*Qf*(x(:,N));
      
      
      % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE 
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      
      ctrl_opt = optimizer(con, obj, sdpsettings('solver','gurobi'), ...
        {x(:,1), xs, us}, u(:,1));
    end
    
    
    % Design a YALMIP optimizer object that takes a position reference
    % and returns a feasible steady-state state and input (xs, us)
    function target_opt = setup_steady_state_target(mpc)

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % INPUTS
      %   ref    - reference to track
      % OUTPUTS
      %   xs, us - steady-state target
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      % Steady-state targets
      n = size(mpc.A,1);
      xs = sdpvar(n, 1);
      us = sdpvar;
      
      % Reference position (Ignore this before Todo 3.2)
      ref = sdpvar;            
            
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE       
      % You can use the matrices mpc.A, mpc.B, mpc.C and mpc.D
      con = [];
      obj = 0;
      % SET PARAMETERS HERE

            
      % CONSTRAINTS
      
      
      % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE 
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % Compute the steady-state target
      target_opt = optimizer(con, obj, sdpsettings('solver', 'gurobi'), ref, {xs, us});
      
    end
  end
end
