close all; clear;

load('CircuitOfAmerica.mat');
load('VehicleDataAudiTTS.mat');
path_ori = preprocessTrack(Track);

%% Extract path information
ds = 5;
path.s = 0:ds:5520;
path.K = path_ori.func_dtheta(path.s);
path.psi = path_ori.func_theta(path.s);
path.cline = path_ori.center(path.s);
path.w_r = path_ori.func_wr(path.s) - 0.8; % shrink the boudary
path.w_l = path_ori.func_wl(path.s) + 0.8; % shrink the boudary

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