% Root finding function via Newton's method
% INPUTS:
%   fun: the function we are computing the root of
%        fun(x) should output [f,dfdx]
%   x0: initial guess for Newton's method
%   dxtol: termination threshold
%   ftol: termination threshold
%   max_iter: maximum iteration limit
%   dxmax: maximum allowed change in x
%
% OUTPUTS:
%   x: estimate for root
%   exit_flag: indicates whether solver succeeded

function [x, exit_flag] = newton_solver(fun,x0,dxtol,ftol,max_iter,dxmax)

    % Initialize exit flag
    exit_flag = 0;

    % Initial root guess
    x_i = x0;

    for n = 1:max_iter

        % Get function value and derivative
        [f_i,df_i] = fun(x_i);

        % Check for zero derivative
        if df_i == 0
            x = x_i;
            exit_flag = -1;
            return;
        end

        % Calculate next Newton guess
        x = x_i - f_i/df_i;

        % Check for excessively large change
        if abs(x - x_i) > dxmax
            exit_flag = -2;
            return;
        end

        % Check change in x
        if abs(x - x_i) <= dxtol
            exit_flag = 1;
            return;
        end

        % Check function value at new guess
        if abs(fun(x)) <= ftol
            exit_flag = 1;
            return;
        end

        % Update current guess
        x_i = x;

    end

end