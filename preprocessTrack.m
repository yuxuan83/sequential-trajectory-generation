function Path = preprocessTrack(Track)
    % options
    VISUALIZE = true;
    CAR_MAKER = true;
    
    % path discretization parameters
    s_window = 50; % moving average window
    ds = 1;
    
    fun_br = @(s) [interp1(Track.arc_s, Track.br(1,:), s, 'spline');
                   interp1(Track.arc_s, Track.br(2,:), s, 'spline');
                   interp1(Track.arc_s, Track.br(3,:), s, 'spline')];

    fun_bl = @(s) [interp1(Track.arc_s, Track.bl(1,:), s, 'spline');
                   interp1(Track.arc_s, Track.bl(2,:), s, 'spline');
                   interp1(Track.arc_s, Track.bl(3,:), s, 'spline')];
               
    arc_s = -100:ds:5600;
    br = fun_br(arc_s);
    bl = fun_bl(arc_s);
    cline = (br + bl) / 2;
    theta = Track.ftheta(arc_s);
    theta(arc_s > 5500) = 5.6;
    theta(arc_s < 0) = -0.67;
    theta = movingAverage(theta, int32(s_window/ds));

    wr = zeros(size(arc_s));
    wl = zeros(size(arc_s));
    for i = 1:length(arc_s)
        wr(i) = norm(cline(:,i) - br(:,i));
        wl(i) = -norm(cline(:,i) - bl(:,i));
    end

    wr = movingAverage(wr, int32(s_window/ds));
    wl = movingAverage(wl, int32(s_window/ds));

    br_avg = zeros(size(br));
    bl_avg = zeros(size(bl));
    for i = 1:length(br_avg)
        br_avg(1:2,i) = cline(1:2,i) + [-wr(i)*sin(theta(i)); wr(i)*cos(theta(i))];
        bl_avg(1:2,i) = cline(1:2,i) + [-wl(i)*sin(theta(i)); wl(i)*cos(theta(i))];
    end
    
    % update left and right boundaries
    br = br_avg;
    bl = bl_avg;
    clear br_avg bl_avg

    % curvature 
    dtheta = zeros(size(arc_s));
    for i = 2:length(arc_s)-1
        dtheta(i) = (theta(i+1)-theta(i-1))/(arc_s(i+1)-arc_s(i-1));
    end
    dtheta(1) = dtheta(2);
    dtheta(end) = dtheta(end-1);

    
    % save data
    Path.arc_s = arc_s;
    Path.cline = cline;
    Path.theta = theta;
    Path.dtheta = dtheta;
    Path.br = br;
    Path.bl = bl;
    Path.wr = wr;
    Path.wl = wl;
    Path.center = @(s) [interp1(Path.arc_s, Path.cline(1,:), s, 'spline');
                        interp1(Path.arc_s, Path.cline(2,:), s, 'spline');
                        interp1(Path.arc_s, Path.cline(3,:), s, 'spline')];
    Path.func_br = @(s) [interp1(Path.arc_s, Path.br(1,:), s, 'spline');
                         interp1(Path.arc_s, Path.br(2,:), s, 'spline');
                         interp1(Path.arc_s, Path.br(3,:), s, 'spline')];
    Path.func_bl = @(s) [interp1(Path.arc_s, Path.bl(1,:), s, 'spline');
                         interp1(Path.arc_s, Path.bl(2,:), s, 'spline');
                         interp1(Path.arc_s, Path.bl(3,:), s, 'spline')];
    Path.func_wr = @(s) interp1(Path.arc_s, Path.wr, s, 'spline');
    Path.func_wl = @(s) interp1(Path.arc_s, Path.wl, s, 'spline');
    Path.func_theta = @(s) interp1(Path.arc, Path.theta, s, 'spline');
    Path.func_dtheta = @(s) interp1(Path.arc, Path.dtheta, s, 'spline');
    
    % save to CarMaker road profile text file
    if CAR_MAKER == true
        arc_s_car_maker = Track.arc_s(11:end-10);
        bl_car_maker = Path.func_br(arc_s_car_maker);
        w_car_maker = Path.func_wr(arc_s_car_maker) - Path.func_wl(arc_s_car_maker);
        
        writematrix(bl_car_maker(1:2,:)', 'track_left_boundary.txt', 'Delimiter', ' ');
        fprintf('track_left_boundary.txt created.\n')
        
        
        % arc length for left boundary
        arc_s_bl_car_maker = zeros(size(arc_s_car_maker));
        arc_s_bl_car_maker(1) = arc_s_car_maker(1);
        
        for i = 2:length(arc_s_bl_car_maker)
            arc_s_bl_car_maker(i) = arc_s_bl_car_maker(i-1) + norm(bl_car_maker(1:2,i)-bl_car_maker(1:2,i-1));
        end
         
        param = ones(length(arc_s_bl_car_maker),1) * [590 -1 0 0 1 0 -999 -999 -999];

        param(:,3) = arc_s_bl_car_maker';
        param(:,6) = w_car_maker';
        
        writematrix(param, 'track_parameters.txt', 'Delimiter', ' ');
        fprintf('track_parameters.txt created.\n')
    end
    
    % plot track
    if VISUALIZE == true
        figure
        plot(br(1,:), br(2,:), 'r.', ...
             bl(1,:), bl(2,:), 'r.', ...
             cline(1,:), cline(2,:), 'b.', ...
             'MarkerSize', 0.1)
        axis equal
        xlabel('$X$ [m]', 'interpreter', 'latex')
        ylabel('$Y$ [m]', 'interpreter', 'latex')

        figure
        subplot(3,1,1)
        plot(arc_s, wr, 'b', ...
             arc_s, wl, 'b', ...
             'LineWidth', 1.5)
        xlim([-inf inf])
        xlabel('$s$ [m]', 'interpreter', 'latex')
        ylabel('$w$ [m]', 'interpreter', 'latex')
        subplot(3,1,2)
        plot(arc_s, Track.ftheta(arc_s), 'k', ...
             arc_s, theta, 'r', ...
             'LineWidth', 1.5)
        xlim([-inf inf])
        xlabel('$s$ [m]', 'interpreter', 'latex')
        ylabel('$\theta$ [m]', 'interpreter', 'latex')
        subplot(3,1,3)
        plot(arc_s, Track.fdtheta(arc_s), 'k', ...
             arc_s, dtheta, 'r', ...
             'LineWidth', 1.5)
        xlim([-inf inf])
        xlabel('$s$ [m]', 'interpreter', 'latex')
        ylabel('$K$ [m]', 'interpreter', 'latex')
    end
end
