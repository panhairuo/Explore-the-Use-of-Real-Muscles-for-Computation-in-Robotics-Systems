 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  constructing network
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% init data structure with default values
data = init_ms_sys_data;  


% change values to fit the task
	
data.num = 20; 			% number of masses
data.show_steps = 1000; % show simulation step every 5000 steps


% range for randomly initialized input weights
data.in_range = [-1 1];  

% defined area for the masses to be places
data.px_lim = [0 10];
data.py_lim = [0 10];

%  defining parameter ranges for
%  spring properties (
%  i.e., nonlinear stiffness and damping functions
data.k_lim = [1 100; 10 100];  
data.d_lim = [1 100; 10 100];



data.show_plot = 0;
data.readout_type = 'LENGTHS'; % using lengths as readout





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loading data for learning and testing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('data_NARMA/NARMA-L2.mat');


% defines data length of learning data, testing data and washout time
wash_out = 100000;		
start = 60000; 
len = 300000;
len_test = 30000;
U = dat.un(start:len,1);  
Y = dat.yn(start:len,1); 

% testing data is taken from another part of the data set
dat_test = dat;  
U_test = dat_test.un(len+1:len+len_test,1);  % use later data for testing
Y_test = dat_test.yn(len+1:len+len_test,1);  


idx=0;  
for i=1:10
 	idx=idx+1;
    disp([' i = ' , num2str(idx)]); 
    %  randomly initialize a network 
    %  with the given parameter
    net = init_ms_sys_net(data); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  simulating network 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic;
    [net2_1,sim_data_1] = simulate_ms_sys_1(net,U);

    % calculate optimal output weights with linear regression
    if (strcmp(net.readout_type,'LENGTHS'))
	    X = sim_data_1.D(wash_out:end,:);  % throw washout away
    else
	    X = sim_data_1.Sx(wash_out:end,:);  % throw washout away
    end

    Yw = Y(wash_out:end,:);
    W_out=X\Yw;


    % start testing with the state right after the learning phase
    net_test_1 = net2_1;

    % set output weights to the optimal ones.
    net_test_1.W_out = W_out;

    % this can be used to check how good the learned, optimal weights
    % represent the learned data
    %  o = X*W_out;
    %  figure;plot(o); hold on;plot(Yw,'r')


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % testing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [net_test_out_1,sim_data_test_1] = simulate_ms_sys_1(net_test_1,U_test);
    t1=toc;


    disp(['MSE_1: ',num2str(mean_squared_error(Y_test,sim_data_test_1.O))])
    fprintf('t1:%.2fs\n',t1)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  simulating network 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic;
    [net2_2,sim_data_2] = simulate_ms_sys_2(net,U);

    % calculate optimal output weights with linear regression
    if (strcmp(net.readout_type,'LENGTHS'))
	    X = sim_data_2.D(wash_out:end,:);  % throw washout away
    else
	    X = sim_data_2.Sx(wash_out:end,:);  % throw washout away
    end

    Yw = Y(wash_out:end,:);
    W_out=X\Yw;


    % start testing with the state right after the learning phase
    net_test_2 = net2_2;

    % set output weights to the optimal ones.
    net_test_2.W_out = W_out;

    % this can be used to check how good the learned, optimal weights
    % represent the learned data
    %  o = X*W_out;
    %  figure;plot(o); hold on;plot(Yw,'r')


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % testing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [net_test_out_2,sim_data_test_2] = simulate_ms_sys_2(net_test_2,U_test);
    t2=toc;


    disp(['MSE_2: ',num2str(mean_squared_error(Y_test,sim_data_test_2.O))])
    fprintf('t2:%.2fs\n',t2)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  simulating network 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic;
    [net2_3,sim_data_3] = simulate_ms_sys_3(net,U);

    % calculate optimal output weights with linear regression
    if (strcmp(net.readout_type,'LENGTHS'))
	    X = sim_data_3.D(wash_out:end,:);  % throw washout away
    else
	    X = sim_data_3.Sx(wash_out:end,:);  % throw washout away
    end

    Yw = Y(wash_out:end,:);
    W_out=X\Yw;


    % start testing with the state right after the learning phase
    net_test_3 = net2_3;

    % set output weights to the optimal ones.
    net_test_3.W_out = W_out;

    % this can be used to check how good the learned, optimal weights
    % represent the learned data
    %  o = X*W_out;
    %  figure;plot(o); hold on;plot(Yw,'r')


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % testing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [net_test_out_3,sim_data_test_3] = simulate_ms_sys_3(net_test_3,U_test);
    t3=toc;


    disp(['MSE_3: ',num2str(mean_squared_error(Y_test,sim_data_test_3.O))])
    fprintf('t3:%.2fs\n',t3)
end