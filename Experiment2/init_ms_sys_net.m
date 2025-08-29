function net = init_ms_sys_net(d)
%  
%   Initializes a mass-spring net with values from a data structure
%  
%   input:	d [OPTIONAL] 	data structure which defines values of the net 
%  							(see init_ms_sys_data in same directory)
%   output: net				ready to simulate mass-spring net
%  
% helmut.hauser@bristol.ac.uk
% For more info see papers:
%  
% Hauser, H.; Ijspeert, A.; Fuchslin, R.; Pfeifer, R. & Maass, W.
% "The role of feedback in morphological computation with compliant bodies"
% Biological Cybernetics, Springer Berlin / Heidelberg, 2012, 106, 595-613
%
% and
%
% Hauser, H.; Ijspeert, A.; Fuchslin, R.; Pfeifer, R. & Maass, W.
% "Towards a theoretical foundation for morphological computation with compliant bodies"
% Biological Cybernetics, Springer Berlin / Heidelberg, 2011, 105, 355-370 
%

if nargin==0
    % With no specific data structure, default values are used
	d = init_ms_sys_data();
end

% Copy init data information to the net structure
net.init_data = d; 

% STRUCTURE of P (information of the point)
net.P.states(:,1) = rand_in_range(d.p_xlim,d.num);
net.P.states(:,2) = rand_in_range(d.p_ylim,d.num);
net.P.states(:,3) = zeros(d.num,1);
net.P.states(:,4) = zeros(d.num,1);
net.P.force       = zeros(d.num,2);  % force acting on the point fx,fy 
net.P.fixed       = zeros(d.num,1);


% Saving initial position and velocities of P (for resetting)
net.init_data.P = net.P;

% Different noises
net.pos_noise = d.pos_noise;  	% noise on the position sensors
net.dist_noise = d.dist_noise;

% Fixing outermost points (convex hull)
k = convhull(net.P.states(:,1), net.P.states(:,2));
net.P.fixed(k,1) = 1; % Fix all points on the convex hull

% Triangulation with Delaunay
tri = delaunay(net.P.states(:,1), net.P.states(:,2));
net.init_data.tri = tri;

% Making a list of connections
tri_s = sort(tri,2); % sorting
R = zeros(d.num,d.num);
for i=1:size(tri,1)
  from = tri_s(i,1); to = tri_s(i,2);  % 1 --> 2
  R(from,to)=1;
  from = tri_s(i,2); to = tri_s(i,3);  % 2 --> 3
  R(from,to)=1;
  from = tri_s(i,1); to = tri_s(i,3);  % 1 --> 3
  R(from,to)=1;
end

w_num = sum(sum(R));
net.W.from = zeros(w_num,1);	  % index of from point
net.W.to   = zeros(w_num,1);    % index of to point
net.W.k1   = rand_in_range_exp(d.k_lim(2,:),w_num); % spring constants
net.W.k3   = rand_in_range_exp(d.k_lim(1,:),w_num); % spring constants
net.W.d1   = rand_in_range_exp(d.d_lim(2,:),w_num); % damping constants
net.W.d3   = rand_in_range_exp(d.d_lim(1,:),w_num); % damping constants
net.W.l0   = zeros(w_num,1);
net.W.dist_old = zeros(w_num,1);
% 在初始化弹簧连接的部分添加
net.W.h = zeros(w_num,1);      % 迟滞变量初始化
net.W.bw_A = 0.1 * ones(w_num,1);    % A参数
net.W.bw_beta = 0.01 * ones(w_num,1); % beta参数
net.W.bw_gamma = 0.01 * ones(w_num,1);% gamma参数
net.W.bw_n = 0.1 * ones(w_num,1);    % n参数
net.W.bw_alpha = 0.05 * ones(w_num,1);% alpha参数

w_idx=0;
for i=1:d.num
   for j=1:d.num
     if(R(i,j)==1)
        w_idx = w_idx+1;
        net.W.from(w_idx,1) = i; 
        net.W.to(w_idx,1) = j;
        net.W.l0(w_idx,1) = e_distance(net.P.states(i,1:2)',net.P.states(j,1:2)');
     end   
   end
end
net.W.dist_old(:,1) = net.W.l0(:,1);

% Input connections
net.W_in = zeros(d.num,d.nInputs);
nmax = ceil(d.num*d.in_conn);
for nI = 1:d.nInputs
   Idx = randperm(d.num);	% random input connections
   Idx(nmax+1:end) = [];
   net.W_in(Idx,nI) = rand_in_range(d.w_in_range,nmax);
end

% Don't allow inputs to fixed points
for i=1:d.num
     if net.P.fixed(i,1)==1
      net.W_in(i,1)=0;
     end
end

% Feedback connections
net.W_fb = zeros(d.num,d.nOutputs);
nmax = ceil(d.num*d.fb_conn);
for nI = 1:d.nOutputs
   Idx = randperm(d.num);
   Idx(nmax+1:end) = [];
   net.W_fb(Idx,nI) = rand_in_range(d.w_fb_range,nmax);	
end

% Output connectivity
switch d.readout_type
    case 'POSITIONS'
        net.W_out = zeros(d.num,d.nOutputs);
        nmax = ceil(d.num*d.out_conn);
        for nO = 1:d.nOutputs
            Odx = randperm(d.num);
            Odx(nmax+1:end) = [];
            net.W_out(Odx,nO) = rand_in_range(d.w_out_range,nmax);
        end
        net.readout_type = d.readout_type;
        
    case 'LENGTHS'
        net.W_out = zeros(size(net.W.l0,1),d.nOutputs);
        nmax = ceil(size(net.W.l0,1)*d.out_conn);
        for nO = 1:d.nOutputs
            Odx = randperm(size(net.W.l0,1));
            Odx(nmax+1:end) = [];
            net.W_out(Odx,nO) = rand_in_range(d.w_out_range,nmax);
        end
        net.readout_type = d.readout_type;
        
    otherwise
        error('Unknown output type chosen');
end

net.fixed_idx = find(net.P.fixed==1);   % indices of fixed points
net.input_idx = find(sum(net.W_in,2)~=0); % indices of input points
net.output_idx = find(sum(net.W_out,2)~=0); % indices of output points

% Show plot if desired
if (d.show_plot == 1)
	plot_graph(net);
end