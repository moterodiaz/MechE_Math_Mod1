%Root finding function via Newton's method
%INPUTS:
%   fun: the function we are computing the root of
%   Note that fun(x) should output [f,dfdx], where dfdx is the derivative of f
%   (see test_func01 below for example)
%   x0: initial guess for Newton's method
%   dxtol: termination threshold (stop when interval abs(x_{i+1}-x_i) < dxtol)
%   ftol: termination threshold (stop when abs(f(x_{i}))<ftol
%   max_iter: maximum iteration limit
%   dxmax: threshold for checking for a divide by zero error: 
%   terminate when abs(x_{i+1}-x_i) > dxmax, where dxmax is a very large number
%OUTPUTS
%   x: estimate for root of fun
%   exit_flag: an integer indicating whether or not the solver succeeded
function [x, exit_flag] = newton_solver(fun,x0,dxtol,ftol,max_iter,dxmax)

    % Initialize exit flag
    exit_flag = 0;

    % Initial root guess
    x_i = x0;

    for n = 1 : max_iter

        % Get function and derivative
        [f_i,df_i] = fun(x_i);

        % Find the next guess
        x = x_i - f_i/df_i;

        % Check for excessively large change
        if abs(x - x_i) > dxmax
            return;
        end
        
        % Check how much the root guess changes
        if abs(x - x_i) <= dxtol
            exit_flag = 1;
            return;
        end
        
        % Check how far the y-value at x is from 0
        if abs(fun(x)) <= ftol
            exit_flag = 1;
            return;
        end

        % If termination criteria isn't met, update x_i
        x_i = x;
    end
    
end