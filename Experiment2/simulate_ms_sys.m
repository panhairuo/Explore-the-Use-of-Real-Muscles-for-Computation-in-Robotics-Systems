function [net_after,sim_data] = simulate_ms_sys (net,input,output)

TEACHER_FORCING = 0; % no teacher forcing is applied 
if nargin == 3
  disp('teacher forcing')
  TEACHER_FORCING = 1;
end

if(isfield(net, 'tansig'))
 disp('using TANSIG springs instead of 3rd order polynom')
end

if(isfield(net, 'rk_steps'))
 disp('using 4th order Runge-Kutta algorithm')
end

% check if we have the case of a symmetric net_after
SYMMETRIC = 0;
THREE_REGIONS = 0;
if (isfield(net,'info')) % right now assuming that symmetric is the only special case we deal here with
  if (strcmp(net.info,'symmetric_net')) 
  	SYMMETRIC = 1;
  	disp('symmetric net !')
  end
  if (strcmp(net.info,'three_regions'))
  	THREE_REGIONS = 1;
  	disp('three regions net !')  
  end
end

% indices of input nodes
in_idx = find(net.W_in~=0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getting necessary out net data structure
P = net.P;
W = net.W;
time_step = net.init_data.time_step;
show_steps = net.init_data.show_steps;
sim_time = size(input,1)*time_step;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num = net.init_data.num;  % for the size of the data matrices
len = size(input,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% --- STABILITY PATCH --- %%%
% Initialize missing states to prevent NaN
if ~isfield(W,'dist_old') || isempty(W.dist_old)
    W.dist_old = zeros(size(W.k1,1),1);
end
if ~isfield(W,'h') || isempty(W.h)
    W.h = zeros(size(W.k1,1),1);
end
% Optional clipping parameters (disabled if not set)
if ~isfield(W,'v_clip');   W.v_clip = Inf;  end     % relative velocity clipping
if ~isfield(W,'h_lim');    W.h_lim  = Inf;  end     % hysteresis variable clipping
eps_h = 1e-12; % prevent |h|^(n-1) singularity at zero
%%% --- STABILITY PATCH END ---

if (net.init_data.save_sim_data==1)
	sim_data.Fx = zeros(len,num);
	sim_data.Fy = zeros(len,num);

	sim_data.Sx_off = zeros(len,num);  % minus the offset
	sim_data.Sy = zeros(len,num);
	sim_data.Sxd = zeros(len,num);
	sim_data.Syd = zeros(len,num);

	sim_data.Sx(1,:)  = P.states(:,1)';	% positions
	sim_data.Sy(1,:)  = P.states(:,2)';
	sim_data.Sxd(1,:) = P.states(:,3)';	% velocities
	sim_data.Syd(1,:) = P.states(:,4)';
end

sim_data.O = zeros(len,net.init_data.nOutputs);
sim_data.D = zeros(len,size(W.k1,1));
sim_data.Sx = zeros(len,num);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Simulation loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idx=0;  
for i=1:len
 	idx=idx+1; 	

 	% set all old forces to zero (to get no unwanted acculumation)
 	P.force(:,1:2) = zeros(num,2);
 	
 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   	% go trough all connections and calculate force
   	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for c=1:size(W.k1,1)
       from = W.from(c,1);
       to   = W.to(c,1);
       p_from = [ P.states(from,1) , P.states(from,2) ]';
       p_to   = [ P.states(to,1) , P.states(to,2) ]';
       [d,ndir] = e_distance( p_from,p_to);

	   if (net.init_data.save_sim_data==1 || strcmp(net.readout_type,'LENGTHS'))
			sim_data.D(idx,c)= d + net.dist_noise*rand(1,1);
	   end

	   if(isfield(net, 'tansig'))
		% ------------ Original tansig branch (unchanged) ------------
	   	A_k = W.k1(c,1);
	   	k_k = W.k3(c,1);
	   	A_d = W.d1(c,1);
	   	k_d = W.d3(c,1);
	   	if(net.tansig == 1)
	   		Fk = -1+(2*A_k)/(1+exp(-2*abs(k_k)*(d-W.l0(c,1)))) + (d-W.l0(c,1));
	   		Fd = -1+(2*A_d)/(1+exp(-2*abs(k_d)*(d-W.dist_old(c,1))/time_step)) + (d-W.dist_old(c,1))/time_step;
	   	else
			Fk = -1+(2*A_k)/(1+exp(-2*abs(k_k)*(d-W.l0(c,1))));
	   		Fd = -1+(2*A_d)/(1+exp(-2*abs(k_d)*(d-W.dist_old(c,1))/time_step));
		end
	  else
		% ------------ Bouc-Wen branch (stability patch) ------------
        %%% --- STABILITY PATCH --- %%%
        % Relative displacement and velocity
        delta_l   = d - W.l0(c,1);                       % only for elastic term
        v_rel     = (d - W.dist_old(c,1)) / time_step;   % for hysteresis/damping

        % Velocity clipping (optional)
        if isfinite(W.v_clip)
            v_rel = max(min(v_rel, W.v_clip), -W.v_clip);
        end

        % Bouc–Wen parameters
        A     = W.bw_A(c,1);
        beta  = W.bw_beta(c,1);
        gamma = W.bw_gamma(c,1);
        n     = W.bw_n(c,1);
        alpha = W.bw_alpha(c,1);
        h     = W.h(c,1);

        % Prevent |h|^(n-1) singularity at h=0,n<1
        abs_h = max(abs(h), eps_h);

        % --- Midpoint method (RK2) integration for h: more stable than explicit Euler ---
        % k1
        dh_dt_1 = A*v_rel - beta*abs(v_rel)*(abs_h^(n-1))*h - gamma*v_rel*(abs_h^n);
        h_half  = h + 0.5*time_step*dh_dt_1;
        abs_hh  = max(abs(h_half), eps_h);
        % k2 (using same v_rel)
        dh_dt_2 = A*v_rel - beta*abs(v_rel)*(abs_hh^(n-1))*h_half - gamma*v_rel*(abs_hh^n);
        h_new   = h + time_step*dh_dt_2;

        % h clipping (optional)
        if isfinite(W.h_lim)
            h_new = max(min(h_new, W.h_lim), -W.h_lim);
        end
        W.h(c,1) = h_new;

        % Unidirectional elasticity (tension only)
        F_el =  W.k3(c,1)*max(0,delta_l)^3 + W.k1(c,1)*max(0,delta_l);
        F_hy =  alpha * W.h(c,1);

        % Damping (Rayleigh type to cubic term), no unidirectional restriction
        Fd   =  W.d1(c,1)*max(0,v_rel) + W.d3(c,1)*(max(0,v_rel)^3);

        % Total force magnitude (along connection direction)
        Fk = F_el + F_hy;
        %%% --- STABILITY PATCH END ---
      end

       % Distribute force to endpoints (direction determined by ndir)
	   if(P.fixed(to,1)==0)   
	   	P.force(to,1)   = P.force(to,1) + (-1)* (Fk+Fd)*ndir(1,1); % f_x
	   	P.force(to,2)   = P.force(to,2) + (-1)* (Fk+Fd)*ndir(2,1); % f_y
	   end
	   if(P.fixed(from,1)==0)   
		P.force(from,1) = P.force(from,1) + (+1)*(Fk+Fd)*ndir(1,1); % f_x
	    P.force(from,2) = P.force(from,2) + (+1)*(Fk+Fd)*ndir(2,1); % f_y
	   end

	   if (net.init_data.save_sim_data==1)
	   	sim_data.Fx(idx,to)   = P.force(to,1);
	   	sim_data.Fy(idx,to)   = P.force(to,2);
	   	sim_data.Fx(idx,from) = P.force(from,1);
	   	sim_data.Fy(idx,from) = P.force(from,2);
	   end

	   % Update historical length
	   W.dist_old(c,1) = d;
    end 
    	   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % add input signals
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if (TEACHER_FORCING==1)
		if (SYMMETRIC || THREE_REGIONS)
			error(['teaching forcing for ***',net.info,'*** net no implemented yet!!'])
		end
    	P.force(:,1) = P.force(:,1) + net.W_in * input(idx,:)' + net.W_fb*output(idx,:)'; 
	else
		if (idx==1)
			if (SYMMETRIC)
				P.force(:,1) = P.force(:,1) + net.W_in(:,1) * input(idx,1)' + net.W_fb*sim_data.O(idx,:)'; % x
				P.force(:,1) = P.force(:,1) + net.W_in(:,2) * input(idx,2)'; % y
			elseif (THREE_REGIONS)
				P.force(:,1) = P.force(:,1) + net.W_in(:,1) * input(idx,1)' + net.W_fb*sim_data.O(idx,:)'; % x
				P.force(:,1) = P.force(:,1) - net.W_in(:,2) * input(idx,2)'; % y-
			else
				P.force(:,1) = P.force(:,1) + net.W_in * input(idx,:)' + net.W_fb*sim_data.O(idx,:)'; % x
			end
		else
			if (SYMMETRIC)
				P.force(:,1) = P.force(:,1) + net.W_in(:,1) * input(idx,1)' + net.W_fb*sim_data.O(idx-1,:)';
				P.force(:,1) = P.force(:,1) + net.W_in(:,2) * input(idx,2)' ;
            elseif(THREE_REGIONS)
				P.force(:,1) = P.force(:,1) + net.W_in(:,1) * input(idx,1)' + net.W_fb*sim_data.O(idx-1,:)';
				P.force(:,1) = P.force(:,1) - net.W_in(:,2) * input(idx,2)' ;
			else
    			P.force(:,1) = P.force(:,1) + net.W_in * input(idx,:)' + net.W_fb*sim_data.O(idx-1,:)';
    		end
   		end
	end
 
    % Fixed point constraints
    P.states(net.fixed_idx,3:4) = zeros(length(net.fixed_idx),2);
    P.force(net.fixed_idx,1:2)  = zeros(length(net.fixed_idx),2);
  
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % integrate mass points
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	for i_node=1:size(P.states,1)
		if(~isfield(net, 'rk_steps'))
	   		x = ode_simple_ms_sys(time_step,P.states(i_node,1:4),P.force(i_node,1:2));
	   		P.states(i_node,:) = x;
	   	else
	   		[tx,xx,vx] = rk4ode2(@rhs_mass_sys, 0, time_step, P.states(i_node,1),P.states(i_node,3), time_step/net.rk_steps,P.force(i_node,1));
	   		[ty,xy,vy] = rk4ode2(@rhs_mass_sys, 0, time_step, P.states(i_node,2),P.states(i_node,4), time_step/net.rk_steps,P.force(i_node,2));
	   		P.states(i_node,1) = xx(1,end); % x-pos
	   		P.states(i_node,2) = xy(1,end); % y-pos
	   		P.states(i_node,3) = vx(1,end); % x-vel
	   		P.states(i_node,4) = vy(1,end); % y-vel 
	   	end

	   	if (net.init_data.save_sim_data==1)
	   		sim_data.Sx(idx,i_node)  = P.states(i_node,1) + net.pos_noise*rand(size(P.states(i_node,1)));
	   		sim_data.Sy(idx,i_node)  = P.states(i_node,2) + net.pos_noise*rand(size(P.states(i_node,2)));
	   		sim_data.Sxd(idx,i_node) = P.states(i_node,3);
	   		sim_data.Syd(idx,i_node) = P.states(i_node,4);
	  	elseif (strcmp(net.readout_type,'POSITIONS'))
	  		sim_data.Sx(idx,i_node) = P.states(i_node,1)+ net.pos_noise*rand(size(P.states(i_node,1)));
	   	end
	end
	
	% Output
	switch net.readout_type
		case 'POSITIONS'
        	sim_data.O(idx,:) = net.W_out' * P.states(:,1); %  assuming just to use x-dimensions        	
       	case 'LENGTHS'
        	sim_data.O(idx,:) = net.W_out' * sim_data.D(idx,:)'; % lengths as outputs	
    end

	if any(isnan(sim_data.O(idx,:)))
		sim_data.ERROR = 'NaN - ERROR / unstable simulation';
		disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
		disp('NaN - ERROR / unstable simulation');
		disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
		net_after = 0;
   		return
    end
end 	

% updating states for the net to be sent back
net_after = net;
net_after.P = P; % dynamic information about the points
net_after.W = W; % dynamic information about the connections

if(isfield(net, 'tansig'))
 net_after.info = 'used TANSIG springs instead of 3rd order polynom';
end

end