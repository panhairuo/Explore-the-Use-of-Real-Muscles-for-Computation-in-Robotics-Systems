 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  constructing network
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% init data structure with default values
data = init_ms_sys_data;  

% change values to fit the task
data.num = 30; 			% number of masses
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

disp('Initializing and testing network...'); 
% randomly initialize a network with the given parameter
net = init_ms_sys_net(data); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
set(gcf,'Color',[1.0 1.0 1.0]); % 深蓝色背景
ax = gca;
set(ax,'Color',[0.2 0.2 0.3],'XColor',[0.8 0.8 1],'YColor',[0.8 0.8 1],...
       'GridColor',[0.4 0.4 0.6],'MinorGridColor',[0.3 0.3 0.4]);

% 绘制目标输出 (亮黄色)
p_target = plot(Y_test,'Color',[1 0.9 0.2],'LineWidth',2.5); 

% hold on;
% % 模型1输出 (绿色)
% p_model1 = plot(sim_data_test_1.O,'--','Color',[0.2 0.8 0.2],'LineWidth',2); 
% % 模型2输出 (橙色)
% p_model2 = plot(sim_data_test_2.O,':','Color',[1 0.6 0.2],'LineWidth',2.2); 
% % 模型3输出 (蓝色)
% p_model3 = plot(sim_data_test_3.O,'-.','Color',[0.2 0.6 1],'LineWidth',2);

% 增强可视化效果
grid on;
box on;
xlim([0 3000]);
title('NARMA (3000 steps)','Color','k','FontSize',16);
xlabel('step','Color','w','FontSize',14);
ylabel('output','Color','w','FontSize',14);

% 自定义图例
leg = legend([p_target, p_model1, p_model2, p_model3],...
    {'target','SD','RB','BH'},...
    'TextColor','w','Location','best');
set(leg,'Color',[0.3 0.3 0.4]);

