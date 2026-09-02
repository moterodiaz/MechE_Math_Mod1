function [x, exit_flag] = secant_solver(fun, x0, x1, dxtol, ftol, max_iter, dxmax)

    x_first = x0;
    x_second = x1;
    denominator_tol = 1e-14;

    f_first = fun(x_first);
    f_second = fun(x_second);

    % check if guess is a root
    if abs(f_first) < ftol
        x = x_first;
        exit_flag = 1;
        return
    elseif abs(f_second) < ftol
        x = x_second;
        exit_flag = 1;
        return
    end

    %iterations until exit
    for iter = 1:max_iter
        denominator = f_second - f_first;

        % make sure denominator passes tolerance
        if ~isfinite(denominator) || abs(denominator) < denominator_tol
            x = x_second;
            exit_flag = -2;
            return
        end

        %secant method algorithm
        dx = -f_second * (x_second - x_first) / denominator;

        % terminate if dx is too large
        if ~isfinite(dx) || abs(dx) > dxmax
            x = x_second;
            exit_flag = -3;
            return
        end

        %secant update
        x_new = x_second + dx;
        f_new = fun(x_new);

        % early termination if near solution
        if abs(dx) < dxtol || abs(f_new) < ftol
            x = x_new;
            exit_flag = 1;
            return
        end

        % shift secant points forward one iteration on update
        x_first = x_second;
        f_first = f_second;

        x_second = x_new;
        f_second = f_new;
    end

    x = x_second;
    exit_flag = 0;
end