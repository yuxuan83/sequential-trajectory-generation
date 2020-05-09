close all; clear; clc;

% Load track and vehicle data
load('CircuitOfAmerica.mat');
load('VehicleData.mat');

%% Extract path information
global mu g
mu = 1.0;
g = 9.81;

s_window = 60; % moving average window
ds = 5;
s = 0:ds:Track.arc_s(end-1);
K = movingAverage(Track.fdtheta(s), int32(s_window/ds));

e_r = Track.center(s) - Track.fun_br(s);
e_l = Track.center(s) - Track.fun_bl(s);
w_r = zeros(size(s));
w_l = zeros(size(s));
for i = 1:length(s)
    w_r(i) = norm(e_r(1:2,i));
    w_l(i) = -norm(e_l(1:2,i));
end

w_r = movingAverage(w_r, int32(s_window/ds));
w_l = movingAverage(w_l, int32(s_window/ds));

path.s = s;
path.K = K;
path.w_r = w_r;
path.w_l = w_l;

figure(1)
subplot(2,1,1)
plot(path.s, path.K, 'k', 'LineWidth', 1.5)
grid on;
xlim([-inf, inf]); 
xlabel('$s$ [m]'); ylabel('$K$ [1/m]')
legend('curvature')
subplot(2,1,2)
plot(path.s, 0*s, '--', ...
     path.s, path.w_r, 'k', ...
     path.s, path.w_l, 'k', ...
     'LineWidth', 1.5)
grid on;
xlim([-inf inf]); ylim([-20 20]); 
xlabel('$s$ [m]'); ylabel('$w$ [m]')
legend('center line', 'right boundary', 'left boundary')

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

% Update speed profile
U_x = calculateSpeedProfile(path, car);

figure(3)
plot(path.s, U_x, 'k', 'LineWidth', 1.5);
grid on;
xlim([-inf inf]); ylim([0, 80]);
xlabel('$s$ [m]'); ylabel('$v$ [m/s]');
legend('speed profile')

% Calculate total lap time
t_s = calculateLapTime(path, U_x);

% Update path
N = length(path.s);
for i = 1:N-1
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
    
    d_k = Ts*[0; ...
            d2; ...
            d3; ...
            d4; ...
            0];
    C = [1 0 0 0 0];
    D = 0;
    [A_k, B_k, ~, ~] = c2dm(A, B, C, D, Ts);
end

%% Function Definitions
function U_x = calculateSpeedProfile(path, car)
    global mu g
    U_x_max = 75;
    
    % 1st pass
    U_x_1 = sqrt(0.9*mu*g./abs(path.K));
    for i = 1:length(U_x_1)
        U_x_1(i) = min(U_x_max, U_x_1(i));
    end
    % 2nd pass
    U_x_2 = U_x_1;
    for i = 1:length(U_x_2)-1       
        ds = path.s(i+1) - path.s(i);
        U_x_2(i+1) = sqrt(U_x_2(i)^2 + 2*(car.F_max/car.m)*ds);
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
        t_s(i) = t_s(i-1) + (path.s(i)-path.s(i-1))/(U_x(i)+U_x(i-1))*2;  
    end
end