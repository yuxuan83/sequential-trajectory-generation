close all; clear; clc

% Load track and vehicle data
load('CircuitOfAmerica.mat');
load('VehicleDataAudiTTS.mat');

global mu g
mu = 0.95;
g = 9.81;

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

figure(1)
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