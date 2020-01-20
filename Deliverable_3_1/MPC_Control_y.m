classdef MPC_Control_y < MPC_Control
  
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
      N = 15;

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
         Q = [0 0 0 0;
           0 0 0 0;
           0 0 0 0;
           0 0 0 1];
          R = 1*eye(m); %50
      end
      
      if (tmp == 1)
          Q = [0 0 0 0;
           0 0 0 0;
           0 0 0 0;
           0 0 0 1];
          R = 50*eye(m); %50
      end
      
      % CONSTRAINTS
      %  u in U = { u | Mu <= m }
      M = [1;-1]; m = [0.3; 0.3];
      % x in X = { x | Fx <= f }
      F = [0 1 0 0; 0 -1 0 0]; f = [0.035; 0.035];
      
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
      
      subplot(321)
      Xf.projection([1 2]).plot();
      xlabel('Angular speed beta [rad/s]');
      ylabel('Angle beta [rad]');
      title('Projection of Terminal invariant set: beta\_dot vs. beta')
      
      subplot(322)
      Xf.projection([1 3]).plot();
      xlabel('Angular speed beta [rad/s]');
      ylabel('Speed x [m/s]');
      title('Projection of Terminal invariant set: beta\_dot vs. x\_dot')
      
      subplot(323)
      Xf.projection([1 4]).plot();
      xlabel('Angular speed beta [rad/s]');
      ylabel('Position x [m]');
      title('Projection of Terminal invariant set: beta\_dot vs. x')
      
      subplot(324)
      Xf.projection([2 3]).plot();
      xlabel('Angle beta [rad]');
      ylabel('Speed x [m/s]');
      title('Projection of Terminal invariant set: beta vs. x\_dot')
      
      subplot(325)
      Xf.projection([2 4]).plot();
      xlabel('Angle beta [rad]');
      ylabel('Position x [m]');
      title('Projection of Terminal invariant set: beta vs. x')
      
      subplot(326)
      Xf.projection([3 4]).plot();
      xlabel('Speed x [m/s]');
      ylabel('Position x [m]');
      title('Projection of Terminal invariant set: x\_dot vs. x')
      end
      
      con = (x(:,2) == mpc.A*(x(:,1)) + mpc.B*(u(:,1))) + (M*u(:,1) <= m);
      obj = (u(:,1))'*R*(u(:,1));
      for i = 2:N-1
          con = con + (x(:,i+1) == mpc.A*(x(:,i)) + mpc.B*(u(:,i)));
          con = con + (F*x(:,i) <= f) + (M*u(:,i) <= m);
          obj = obj + (x(:,i))'*Q*(x(:,i)) + (u(:,i))'*R*(u(:,i));
      end
      con = con + (Ff*x(:,N) <= ff);
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
