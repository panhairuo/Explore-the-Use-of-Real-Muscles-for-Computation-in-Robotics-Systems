clear;
%% 1. Data Preparation - Generate robot trajectory and torque
% Robot parameters
l1 = 0.5; l2 = 0.5; % arm length (m)
m1 = 0.1; m2 = 0.1; % mass (kg)
g = 0;           % gravity acceleration

%% ---------------- Trajectory: First quadrant figure-8 ----------------
t = 0:0.001:40; % 40-second trajectory, 1000Hz sampling
dt = t(2) - t(1);

A = 0.4;  % x amplitude
B = 0.2;  % y amplitude
freq = 0.25; % base frequency

% Figure-8 trajectory (Lissajous variant)
x = A * sin(2*pi*freq*t);
y = B * sin(4*pi*freq*t); % y frequency is twice x

% Translate to first quadrant (keep 0.1 m safety margin)
x = x - min(x) + 0.1;
y = y - min(y) + 0.3;

% Scale to avoid exceeding arm reach
max_reach = 0.95 * (l1 + l2);
r = sqrt(x.^2 + y.^2);
if any(r > max_reach)
    scale_factor = max_reach / max(r);
    x = x * scale_factor;
    y = y * scale_factor;
end

figure('Name','Trajectory');
plot(x, y, 'LineWidth', 1.5);
title('First Quadrant Figure-8 Trajectory');
xlabel('X Position (m)');
ylabel('Y Position (m)');
axis equal; grid on;

% ===== 1. Inverse Kinematics =====
theta1 = zeros(size(t));
theta2 = zeros(size(t));
dtheta1 = zeros(size(t));
dtheta2 = zeros(size(t));
ddtheta1 = zeros(size(t));
ddtheta2 = zeros(size(t));

for i = 1:length(t)
    r = sqrt(x(i)^2 + y(i)^2);
    
    cos_theta2 = (r^2 - l1^2 - l2^2) / (2*l1*l2);
    cos_theta2 = min(max(cos_theta2, -1), 1);  % prevent out-of-bound
    theta2(i) = acos(cos_theta2);              % elbow-down solution
    
    k1 = l1 + l2*cos(theta2(i));
    k2 = l2*sin(theta2(i));
    theta1(i) = atan2(y(i), x(i)) - atan2(k2, k1);
end

theta1 = unwrap(theta1);
theta2 = unwrap(theta2);

% Angular velocity
dtheta1(2:end) = diff(theta1) / dt;
dtheta2(2:end) = diff(theta2) / dt;

% Angular acceleration
ddtheta1(2:end) = diff(dtheta1) / dt;
ddtheta2(2:end) = diff(dtheta2) / dt;

% ===== 2. Inverse Dynamics =====
tau1 = zeros(size(t));
tau2 = zeros(size(t));

for i = 1:length(t)
    th1 = theta1(i);
    th2 = theta2(i);
    dth1 = dtheta1(i);
    dth2 = dtheta2(i);
    ddth1 = ddtheta1(i);
    ddth2 = ddtheta2(i);
    
    % Inertia matrix M(θ)
    M = [m1*l1^2 + m2*(l1^2 + l2^2 + 2*l1*l2*cos(th2)),   m2*(l2^2 + l1*l2*cos(th2));
         m2*(l2^2 + l1*l2*cos(th2)),                      m2*l2^2];
    
    % Coriolis and centrifugal forces C(θ, θ̇)
    C = [-m2*l1*l2*sin(th2)*(2*dth1*dth2 + dth2^2);
          m2*l1*l2*sin(th2)*(dth1^2)];
    
    % Gravity term G(θ)
    G = [(m1 + m2)*g*l1*cos(th1) + m2*g*l2*cos(th1+th2);
          m2*g*l2*cos(th1+th2)];
    
    % Theoretical inverse dynamics torque
    tau = M * [ddth1; ddth2] + C + G;
    tau1(i) = tau(1);
    tau2(i) = tau(2);
end

% Split into training and testing sets (70% training, 30% testing)
train_ratio = 0.7;
train_len = round(length(t)*train_ratio);
U_train = [x(1:train_len)', y(1:train_len)'];
Y_train = [tau1(1:train_len)', tau2(1:train_len)'];

U_test = [x(train_len+1:end)', y(train_len+1:end)'];
Y_test = [tau1(train_len+1:end)', tau2(train_len+1:end)'];

%% 2. Bouc-Wen Model
data = init_ms_sys_data;
data.num = 30;
data.show_steps = 1000;
data.in_range = [-1 1];
data.px_lim = [0 10];
data.py_lim = [0 10];
data.k_lim = [1 100; 10 100];
data.d_lim = [10 100; 1 10];
data.show_plot = 0;
data.readout_type = 'POSITIONS';
data.nInputs = 2;
data.nOutputs = 2;

net.rk_steps = 4;
net.W.v_clip = 5;
net.W.h_lim  = 1e3;
net = init_ms_sys_net(data);

num_springs = length(net.W.k1);
net.W.bw_A = 1.0 * ones(num_springs,1);
net.W.bw_beta = 0.1 * ones(num_springs,1);
net.W.bw_gamma = 0.1 * ones(num_springs,1);
net.W.bw_n = 1.0 * ones(num_springs,1);
net.W.bw_alpha = 0.5 * ones(num_springs,1);
net.W.h = zeros(num_springs,1);

%% 3. Train Bouc-Wen Model
tic;
[net_train, sim_data_train] = simulate_ms_sys(net, U_train);

wash_out = 8000;
if (strcmp(net.readout_type,'LENGTHS'))
    X_train = sim_data_train.D(wash_out:end,:);
else
    X_train = sim_data_train.Sx(wash_out:end,:);
end
W_out = X_train \ Y_train(wash_out:end,:);
net_test = net_train;
net_test.W_out = W_out;

% Test
[~, sim_data_test] = simulate_ms_sys(net_test, U_test);
bw_pred = sim_data_test.O;
bw_time = toc;

% MSE
bw_mse1 = mean((bw_pred(:,1) - Y_test(:,1)).^2);
bw_mse2 = mean((bw_pred(:,2) - Y_test(:,2)).^2);
fprintf('Bouc-Wen Model MSE - Tau1: %.4e, Tau2: %.4e, Time: %.2fs\n', ...
        bw_mse1, bw_mse2, bw_time);

%% 4. Linear Regression (LR) Model
tic;
U_train_lr = [U_train, ones(size(U_train,1),1)];
U_test_lr  = [U_test,  ones(size(U_test,1),1)];

W_lr = U_train_lr \ Y_train;
lr_pred = U_test_lr * W_lr;
lr_time = toc;

lr_mse1 = mean((lr_pred(:,1) - Y_test(:,1)).^2);
lr_mse2 = mean((lr_pred(:,2) - Y_test(:,2)).^2);
fprintf('LR Model MSE - Tau1: %.4e, Tau2: %.4e, Time: %.4fs\n', ...
        lr_mse1, lr_mse2, lr_time);

%% 5. ANN Model
tic;
inputSize = 2;       % Input features: x,y
numHiddenUnits = 30; % Hidden layer neurons
numResponses = 2;    % Output: tau1, tau2

% Define ANN network structure
layers = [ ...
    featureInputLayer(inputSize)
    fullyConnectedLayer(numHiddenUnits)
    reluLayer
    fullyConnectedLayer(numResponses)
    regressionLayer];

options = trainingOptions('adam', ...
    'MaxEpochs',100, ...       % adjustable
    'MiniBatchSize', 256, ...
    'InitialLearnRate',1e-3, ...
    'Shuffle','every-epoch', ...
    'Verbose',0);

% Training
net_ann = trainNetwork(U_train, Y_train, layers, options);

% Prediction
ann_pred = predict(net_ann, U_test);
ann_time = toc;

% MSE
ann_mse1 = mean((ann_pred(:,1) - Y_test(:,1)).^2);
ann_mse2 = mean((ann_pred(:,2) - Y_test(:,2)).^2);
fprintf('ANN Model MSE - Tau1: %.4e, Tau2: %.4e, Time: %.2fs\n', ...
        ann_mse1, ann_mse2, ann_time);

%% 6. Visualization Comparison
figure;
subplot(2,1,1);
plot(Y_test(:,1), 'r', 'LineWidth', 1.5); hold on;
plot(lr_pred(:,1), 'k:');
plot(bw_pred(:,1), 'b--');
plot(ann_pred(:,1), 'm-.');
title('Joint 1 Torque Prediction Comparison');
legend('True', 'LR', 'Bouc-Wen', 'ANN');
grid on;

subplot(2,1,2);
plot(Y_test(:,2), 'r', 'LineWidth', 1.5); hold on;
plot(lr_pred(:,2), 'k:');
plot(bw_pred(:,2), 'b--');
plot(ann_pred(:,2), 'm-.');
title('Joint 2 Torque Prediction Comparison');
legend('True', 'LR', 'Bouc-Wen', 'ANN');
grid on;

%% 7. MSE Bar Chart
figure;
mse_all = [lr_mse1, bw_mse1, ann_mse1;
           lr_mse2, bw_mse2, ann_mse2];
bar(mse_all);
set(gca, 'XTickLabel', {'Joint 1', 'Joint 2'});
legend('LR', 'Bouc-Wen', 'ANN');
ylabel('MSE');
title('Performance Comparison of Different Models');
grid on;

%% 8. Time Bar Chart
figure;
time_all = [lr_time, bw_time, ann_time];
bar(time_all);
set(gca, 'XTickLabel', {'LR', 'Bouc-Wen', 'ANN'});
ylabel('Time (s)');
title('Runtime Comparison of Different Models');
grid on;