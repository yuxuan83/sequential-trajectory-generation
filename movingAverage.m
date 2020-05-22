function x_avg = movingAverage(x, n)
%%% Moving average 
%%% x: vector, n: window size
    if mod(n,2) ~= 0
        n = n-1;
    end
    
    a = n/2;
    x_avg = zeros(size(x));
    for i = 1+a:length(x)-a
        x_avg(i) = mean(x(i-a:i+a));
    end
    x_avg(1:a) = x_avg(1+a);
    x_avg(end-a+1:end) = x_avg(end-a);
end