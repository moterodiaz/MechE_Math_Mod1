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
function [x, exit_flag] = bisection_solver(fun,x_left,x_right,dxtol,ftol,max_iter) 
 
    % Initialize exit flag 
    exit_flag = 0; 
 
    if sign(fun(x_left)) == sign(fun(x_right)) 
        disp('Root Does Not Necessarily Exist Between Bounds') 
        return; 
    end 
 
    % Loop function max amount of times if necessary 
    for i = 1 : max_iter 
        
    % Calculate Midpoint 
        c = (x_left + x_right)/2; 
        x = c; 
 
    % Calculate function at midpoint 
    f_c = fun(c); 
     
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
    if sign(fun(x_left)) ~= sign(f_c) 
        x_right = c; 
    else 
        x_left = c; 
    end 
     
    % Termination threshold --> Diff between left and right bounds 
        if abs(x_left - x_right) <= dxtol 
            exit_flag = 1; 
            disp('X termination criteria met'); 
            return; 
        end 
 
    end