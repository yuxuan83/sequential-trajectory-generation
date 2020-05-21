close all; clear; clc;

% Load track and vehicle data
load('CircuitOfAmerica.mat');
load('VehicleDataAudiTTS.mat');
path_ori = preprocessTrack(Track);

%% Extract path information
global mu g
mu = 0.95;
g = 9.81;

ds = 2.75;
path.s = 0:ds:5520;
path.K = path_ori.func_dtheta(path.s);
path.psi = path_ori.func_theta(path.s);
path.cline = path_ori.center(path.s);
path.w_r = path_ori.func_wr(path.s);
path.w_l = path_ori.func_wl(path.s);

figure(1)
subplot(2,1,1)
plot(path.s, path.K, 'k', 'LineWidth', 1.5)
grid on;
xlim([-inf, inf]); 
xlabel('$s$ [m]', 'interpreter', 'latex'); ylabel('$K$ [1/m]', 'interpreter', 'latex')
legend('curvature', 'interpreter', 'latex')
subplot(2,1,2)
plot(path.s, 0*path.s, 'b-', ...
     path.s, path.w_r, 'k', ...
     path.s, path.w_l, 'k', ...
     'LineWidth', 1.5)
grid on;
xlim([-inf inf]); ylim([-20 20]); 
xlabel('$s$ [m]', 'interpreter', 'latex'); ylabel('$w$ [m]', 'interpreter', 'latex')
legend('center line', 'right boundary', 'left boundary', 'interpreter', 'latex')

%% Fiala tire model table
dalpha = 0.0001;
diff_num = @(f,x,dx) (f(x+dx)-f(x-dx))/(2*dx);

% front tire
F_z_f = car.m*g*car.b/(car.a+car.b);
alpha_f_cr = atan(3*mu*F_z_f/car.C_f);
alpha_f = -alpha_f_cr:dalpha:alpha_f_cr;
F_y_f = fialaModel(alpha_f, car.C_f, F_z_f, mu);

% rear tire
F_z_r = car.m*g*car.a/(car.a+car.b);
alpha_r_cr = atan(3*mu*F_z_r/car.C_r);
alpha_r = -alpha_r_cr:dalpha:alpha_r_cr;
F_y_r = fialaModel(alpha_r, car.C_r, F_z_r, mu);

% alpha_tilda function
alpha_f_tilda_func = @(F_y) interp1(F_y_f, alpha_f, F_y);
alpha_r_tilda_func = @(F_y) interp1(F_y_r, alpha_r, F_y);

% C_tilda function
F_y_f_func = @(alpha) interp1(alpha_f, F_y_f, alpha, 'pchip', 'extrap');
F_y_r_func = @(alpha) interp1(alpha_r, F_y_r, alpha, 'pchip', 'extrap');
C_f_tilda_func = @(alpha) diff_num(F_y_f_func, alpha, dalpha);
C_r_tilda_func = @(alpha) diff_num(F_y_r_func, alpha, dalpha);

F_y_exp = 4000;
alpha_f_tilda = alpha_f_tilda_func(F_y_exp);
alpha_r_tilda = alpha_r_tilda_func(F_y_exp);
C_f_tilda = C_f_tilda_func(alpha_f_tilda);
C_r_tilda = C_r_tilda_func(alpha_r_tilda);

figure(2)
plot(alpha_f, F_y_f, 'r', ...   
     alpha_r, F_y_r, 'b', ...
     alpha_f, C_f_tilda*(alpha_f-alpha_f_tilda)+F_y_exp, 'g', ...
     alpha_r, C_r_tilda*(alpha_r-alpha_r_tilda)+F_y_exp, 'g', ...
     alpha_f_tilda, F_y_exp, 'g*', ...
     alpha_r_tilda, F_y_exp, 'go', ...
     'LineWidth', 1.5)
grid on;
xlim([-inf inf]); ylim([-1e4, 1e4]);
xlabel('$\alpha$ [rad]'); ylabel('$F_y$ [N]');
title('Lateral Force')
legend('Front', 'Rear')

%% main function loop

%% Update speed profile
tic
fprintf('Speed Profile Update ... ');
U_x = calculateSpeedProfile(path, car);
fprintf(' %.2fs\n', toc)

figure(3)
plot(path.s, U_x, 'k', 'LineWidth', 1.5);
grid on;
xlim([-inf inf]); ylim([0, 80]);
xlabel('$s$ [m]'); ylabel('$v$ [m/s]');
title('speed profile')

%% Calculate total lap time
t_s = calculateLapTime(path, U_x);

%% Update path
tic
fprintf('Path Update ... ');

N = length(path.s);
% Construct cost function
% initialize H_sr
lambda = 1;
H_sr = zeros(2*(N-1), 6*N);
for i = 1:N-1
    % define current row and column indices
    row_idx = 2*(i-1);
    col_idx = 6*(i-1);
    
    ds = path.s(i+1) - path.s(i);
    H_k = [0, 0, 0, 0, 1/ds, 0; ...
           0, 0, 0, 0, 0, sqrt(lambda)];
    
    H_sr(row_idx+1:row_idx+2, col_idx+1:col_idx+6) = -H_k;
    H_sr(row_idx+1:row_idx+2, col_idx+7:col_idx+12) = H_k;
end

% Construct equality constraints 15b
% initialize A_eq and b_eq
A_eq = zeros(5*(N+1), 6*N);
b_eq = zeros(5*(N+1), 1);

% 1nd - N-1 row
for i = 1:N-1
    % construct time varying model of each step
    Ts = t_s(i+1) - t_s(i);
    
    F_y_f_tilda = car.m*car.b/(car.a+car.b)*U_x(i)^2*path.K(i);
    F_y_r_tilda = car.m*car.a/(car.a+car.b)*U_x(i)^2*path.K(i);
    alpha_f_tilda = alpha_f_tilda_func(F_y_f_tilda);
    alpha_r_tilda = alpha_r_tilda_func(F_y_r_tilda);
    C_f_tilda = C_f_tilda_func(alpha_f_tilda);
    C_r_tilda = C_r_tilda_func(alpha_r_tilda);
    
    a33 = -(car.a^2*C_f_tilda + car.b^2*C_r_tilda) / (U_x(i)*car.I);
    a34 = (car.b*C_r_tilda - car.a*C_f_tilda) / car.I;
    a43 = (car.b*C_r_tilda - car.a*C_f_tilda) / (car.m*U_x(i)^2) - 1;
    a44 = -(C_f_tilda + C_r_tilda) / (car.m*U_x(i));
    
    A = [0 U_x(i) 0 U_x(i) 0; ...
         0 0 1 0 0; ...
         0 0 a33 a34 0; ...
         0 0 a43 a44 0; ...
         0 0 1 0 0];
    B = [0; ...
         0; ...
         car.a*C_f_tilda/car.I; ...
         C_f_tilda/(car.m*U_x(i)); ...
         0]; 
    
    d2 = -path.K(i)*U_x(i);
    d3 = (car.a*C_f_tilda*alpha_f_tilda - car.b*C_r_tilda*alpha_r_tilda + car.a*F_y_f_tilda - car.b*F_y_r_tilda) / car.I;
    d4 = (C_f_tilda*alpha_f_tilda + C_r_tilda*alpha_r_tilda + F_y_f_tilda + F_y_r_tilda) / (car.m*U_x(i));
    
    d_k = [ 0; ...
           d2; ...
           d3; ...
           d4; ...
            0]*Ts;
    C = [1 0 0 0 0];
    D = 0;
    [A_k, B_k, ~, ~] = c2dm(A, B, C, D, Ts);
    
    % define current row and column indices
    row_idx = 5*(i-1);
    col_idx = 6*(i-1);
    
    A_eq(row_idx+1:row_idx+5, col_idx+7:col_idx+12) = [eye(5), zeros(5,1)];    % block (i,i+1)
    A_eq(row_idx+1:row_idx+5, col_idx+1:col_idx+6) = -[A_k, B_k];              % block (i,i)
    b_eq(row_idx+1:row_idx+5, 1) = d_k;                                        % block (i,1)
end

% N-th row: initial condition
A_eq(5*(N-1)+1:5*(N-1)+5, 1:6) = [eye(5), zeros(5,1)];                         % block (N,1)
b_eq(5*(N-1)+1:5*(N-1)+5, 1) = [0; 0; 0; 0; path.psi(1)];                      % block (N,1)

% (N+1)-th row: final condition
A_eq(5*N+1:5*N+5, 6*(N-1)+1:6*(N-1)+6) = [diag([1 0 1 0 1]), zeros(5,1)];      % block (N+1,N)
b_eq(5*N+1:5*N+5, 1) = [0; 0; 0; 0; path.psi(end)];                            % block (N+1,1)

% Construct inequality contraints 15c
% initialize A_in and b_in
A_in = zeros(8*N, 6*N);
b_in = zeros(8*N, 1);

delta_max = 37*pi/180;
delta_min = -delta_max;
for i = 1:N
    G_k = [1, 0, 0, 0, 0, 0; ...
           0, 0, car.a/U_x(i), 1, 0,-1; ...
           0, 0,-car.b/U_x(i), 1, 0, 0; ... 
           0, 0, 0, 0, 0, 1];
    b_in_k_max = [path.w_r(i)-1; ...
                  alpha_f_cr; ...
                  alpha_r_cr; ...
                  delta_max];
    b_in_k_min = [ path.w_l(i)+1; ...
                  -alpha_f_cr; ...
                  -alpha_r_cr; ...
                   delta_min];
    
    row_idx = 8*(i-1);
    col_idx = 6*(i-1);
    A_in(row_idx+1:row_idx+8, col_idx+1:col_idx+6) = [ G_k; ...    % block (i,i)
                                                      -G_k];
    b_in(row_idx+1:row_idx+8, 1) = [ b_in_k_max; ...
                                    -b_in_k_min];                  % block (i,1)
end

% Solve the convex optimization problem using CVX
n = size(H_sr,2);
cvx_begin quiet
    cvx_precision low
    variable x_aug_star(n);
    minimize( norm(H_sr*x_aug_star,2) )
    subject to
        A_eq*x_aug_star == b_eq
        A_in*x_aug_star <= b_in
cvx_end
fprintf(' %.2fs\n', toc);

% Retreive the original states
x_aug_star = reshape(x_aug_star, 6, length(x_aug_star)/6);

e_star = x_aug_star(1,:);
delta_psi_star = x_aug_star(2,:);
omega_star = x_aug_star(3,:);
beta_star = x_aug_star(4,:);
psi_star = x_aug_star(5,:);
delta_star = x_aug_star(6,:);

X = path.cline(1,:) - e_star.*sin(path.psi);
Y = path.cline(2,:) + e_star.*cos(path.psi);

%% Plot results
figure(4)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1])
clf; hold on;
plot(X, Y, 'r', 'LineWidth', 1.5)
plot(path_ori.bl(1,:), path_ori.bl(2,:), '-.k', ...
     path_ori.br(1,:), path_ori.br(2,:), '-.k', ...
     'MarkerSize', 0.4, 'LineWidth', 0.4);
plot(path_ori.cline(1,:), path_ori.cline(2,:), '-.k',...
     'MarkerSize', 0.4, 'LineWidth', 0.4)
grid on;
xlim([-400 1500]); ylim([-500 1400]); axis square
xlabel('$X$ [m]', 'interpreter', 'latex'); ylabel('$Y$ [m]', 'interpreter', 'latex');
title('Track', 'interpreter', 'latex')
legend('Optimized Path')

figure(5)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1])
subplot(6,1,1)
plot(path.s, e_star, 'LineWidth', 1.5);
xlim([-inf inf])
ylabel('$e^*$ [m]', 'interpreter', 'latex')
title('Lateral Errors', 'interpreter', 'latex')
subplot(6,1,2)
plot(path.s, delta_psi_star, 'LineWidth', 1.5);
xlim([-inf inf])
ylabel('$\Delta\psi^*$ [rad]')
title('Heading Angle Errors', 'interpreter', 'latex');
subplot(6,1,3)
plot(path.s, omega_star, 'LineWidth', 1.5);
xlim([-inf inf])
ylabel('$\dot{\psi^*}$ [rad/s]')
title('Steering Angle', 'interpreter', 'latex')
subplot(6,1,4)
plot(path.s, beta_star, 'LineWidth', 1.5)
xlim([-inf inf])
ylabel('$\beta^*$ [rad]', 'interpreter', 'latex')
title('Side Slip Angle [rad]', 'interpreter', 'latex')
subplot(6,1,5)
plot(path.s, psi_star, ...
     path.s, path.psi , ...
     'LineWidth', 1.5)
xlim([-inf inf])
xlabel('$s$ [m]', 'interpreter', 'latex'); ylabel('$\psi^*$ [rad]', 'interpreter', 'latex')
title('Heading Angle', 'interpreter', 'latex')
legend('Optimized Heading Angle', 'Original Path Heading Angle', 'interpreter', 'latex')
subplot(6,1,6)
plot(path.s, delta_star, 'LineWidth', 1.5)
xlim([-inf inf])
xlabel('$s$ [m]', 'interpreter', 'latex'); ylabel('$\delta^*$ [rad]', 'interpreter', 'latex')

%% Function Definitions
function U_x = calculateSpeedProfile(path, car)
    global mu g
    U_x_max = car.v_max;
    
    % 1st pass
    U_x_1 = sqrt(0.95*mu*g./abs(path.K));
    for i = 1:length(U_x_1)
        U_x_1(i) = min(U_x_max, U_x_1(i));
    end
    % 2nd pass
    U_x_2 = U_x_1;
    for i = 1:length(U_x_2)-1       
        ds = path.s(i+1) - path.s(i);
        F_max = car.F_max - 0.5*car.C_d*car.rho*car.A*U_x_2(i);
        U_x_2(i+1) = sqrt(U_x_2(i)^2 + 2*(F_max/car.m)*ds);
        U_x_2(i+1) = min(U_x_2(i+1), U_x_1(i+1));
    end
    % 3rd pass
    U_x_3 = U_x_2;
    for i = length(U_x_3)-1:-1:1
        ds = path.s(i+1) - path.s(i);
        U_x_3(i) = sqrt(U_x_3(i+1)^2 - 2*(car.F_min/car.m)*ds);
        U_x_3(i) = min(U_x_3(i), U_x_2(i));
    end

    U_x = U_x_3;
end

function t_s = calculateLapTime(path, U_x) 
    t_s = zeros(size(path.s));
    for i = 2:length(t_s)
        t_s(i) = t_s(i-1) + (path.s(i)-path.s(i-1)) / ((U_x(i)+U_x(i-1))/2);
    end
end