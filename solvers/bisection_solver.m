%Root finding function via bisection algorithm 
%INPUTS: 
%   fun: the function we are computing the root of 
%   x_left: left guess 
%   x_right: right guess 
%   note that f(x_left) and f(x_right) should have different signs 
%   dxtol: termination threshold (stop when interval x_right-x_left < dxtol) 
%   ftol: termination threshold (stop when abs(f(x_guess))<ftol 
%   max_iter: maximum iteration limit 
%OUTPUTS
%   x: estimate for root of fun
%   exit_flag: an integer indicating whether or not the solver succeeded
%       1: converged (ftol or dxtol met)
%       0: max_iter reached
%      -1: no sign change between x_left and x_right
%      -2: non-finite function value at a midpoint
%   guess_list: the endpoint discarded at each iteration (cleaner than the
%   midpoints for convergence plots, since it doesn't get lucky and land
%   near the root early)
function [x, exit_flag, guess_list] = bisection_solver(fun,x_left,x_right,dxtol,ftol,max_iter)

    % Initialize exit flag
    exit_flag = 0;
    guess_list = [];
    x = NaN;

    % Evaluate the endpoints once; f_left is kept current inside the loop
    f_left = fun(x_left);
    f_right = fun(x_right);

    % Check if either bound is already a root
    if abs(f_left) <= ftol
        x = x_left;
        exit_flag = 1;
        return;
    elseif abs(f_right) <= ftol
        x = x_right;
        exit_flag = 1;
        return;
    end

    if sign(f_left) == sign(f_right)
        exit_flag = -1;
        return; 
    end 
 
    % Loop function max amount of times if necessary 
    for i = 1 : max_iter 
        
        % Calculate Midpoint 
        c = (x_left + x_right)/2; 
        x = c; 
 
        % Calculate function at midpoint 
        f_c = fun(c); 

        if ~isfinite(f_c)
            exit_flag = -2;
            return;
        end
     
        % Termination threshold --> Diff between estimate and zero 
        if abs(f_c) <= ftol 
            exit_flag = 1; 
            return; 
        end 
 
        % Create new bounds using midpoint 
         
        % If the left bound and midpoint have different sign values, midpoint 
        % becomes new right bound 
        % If the right bound and midpoint have different sign values, midpoint 
        % becomes new left bound 
        if sign(f_left) ~= sign(f_c)
            guess_list(end+1) = x_right; % x_right is being discarded
            x_right = c;
        else
            guess_list(end+1) = x_left; % x_left is being discarded
            x_left = c;
            f_left = f_c;
        end
     
        % Termination threshold --> Diff between left and right bounds 
        if abs(x_left - x_right) <= dxtol 
            exit_flag = 1; 
            return; 
        end 
 
    end
end
