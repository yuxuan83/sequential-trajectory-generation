close all; clear; clc

load('CircuitOfAmerica.mat')

bl = Track.bl(1:2,11:587)';
writematrix(bl, 'track_lb.txt', 'Delimiter', ' ');

arc_s = 0:5500;
bl = Track.fun_bl(arc_s);
bl = bl(1:2,:)';
br = Track.fun_br(arc_s);
br = br(1:2,:)';

v = br - bl;

s = 0;
arc_s_l = zeros(size(arc_s));
for i = 2:length(v)
    s = s + norm(bl(i,:)-bl(i-1,:));
    arc_s_l(i) = s;
end

for i = 1:length(v)
    w_ori(i,1) = norm(v(i,:));
end

w = movingAverage(w_ori, 60);
   
% export 
param = [];
for i = 1:length(v)
    
    if mod(i,10) == 0
        param = [param; 
                 590 -1 arc_s_l(i) 0 1 w(i) -999 -999 -999];
    end
end

writematrix(param, 'track_param.txt', 'Delimiter', ' ');


