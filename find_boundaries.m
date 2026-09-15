function Q = find_boundaries(x0,y0,theta, egg_params)
    dxtol = 1e-14;
    ftol = 1e-14;
    max_iter = 100;
    dxmax = 10;

    Vver = zeros(2,1);
    Vhor = zeros(2,1);
    x_hor = zeros(2,1);
    y_hor = zeros(2,1);
    x_ver = zeros(2,1);
    y_ver = zeros(2,1);
   
    h1 = @(s) egg_wrapper_x1 (s, x0, y0, theta, egg_params) ;
    h2 = @(s) egg_wrapper_y1 (s, x0, y0, theta, egg_params) ;

    %compute the value of s for which the corresponding point on the oval
    %has an x-coordinate of zero
    s_rootx = secant_solver(h1,0,.1, dxtol, ftol, max_iter, dxmax);
    s_rooty = secant_solver(h2,0,.1, dxtol, ftol, max_iter, dxmax);
    
    [V,~] = egg_func(s_rootx,x0,y0,theta,egg_params);
    x_hor(1) = V(1);
    y_hor(1) = V(2);

    [V,~] = egg_func(s_rooty,x0,y0,theta,egg_params);
    x_ver(1) = V(1);
    y_ver(1) = V(2);

    hor = true;
    ver = true;

    for i=0:0.25:1
        if ~hor && ~ver
            break;
        end

        if hor
            s_rootx = secant_solver(h1,i,i+.1, dxtol, ftol, max_iter, dxmax);
            [Vhor,~] = egg_func(s_rootx,x0,y0,theta,egg_params);

            cond1 = (Vhor(1) ~= x_hor(1)) && (Vhor(2) ~= y_hor(1));
            if cond1
                x_hor(2) = Vhor(1);
                y_hor(2) = Vhor(2);
                hor = false;
            end
        end

        if ver
            s_rooty = secant_solver(h2,i,i+.1, dxtol, ftol, max_iter, dxmax);
            [Vver,~] = egg_func(s_rooty,x0,y0,theta,egg_params);

            cond2 = (Vver(1) ~= x_ver(1)) && (Vver(2) ~= y_ver(1));
            if cond2()
                x_ver(2) = Vver(1);
                y_ver(2) = Vver(2);
                ver = false;
            end
        end
    end

    x_l = min(x_hor); x_r = max(x_hor);
    y_l = min(y_hor); y_r = max(y_hor);
    x_d = min(x_ver); x_u = max(x_ver);
    y_d = min(y_ver); y_u = max(y_ver);
   

    Q = [x_l, y_l; x_r, y_r; x_d, y_d; x_u, y_u];
    array2table(Q, "VariableNames",{'x', 'y'}, "RowNames", {'Left', 'Right', 'Down', 'Up'})
end


%wrapper function that calls egg_func
%and only returns the x coordinate of the
%point on the perimeter of the egg
%(single output)
function x_out = egg_wrapper_x1(s,x0,y0,theta,egg_params)
    [~, G] = egg_func(s,x0,y0,theta,egg_params);
    x_out = G(1);
end
function x_out = egg_wrapper_y1(s,x0,y0,theta,egg_params)
    [~, G] = egg_func(s,x0,y0,theta,egg_params);
    x_out = G(2);
end