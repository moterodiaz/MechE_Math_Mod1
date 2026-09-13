% Ensure all repo subfolders are in MATLAB's path
root_dir = fileparts(mfilename('fullpath'));
addpath(genpath(root_dir));

global input_list;
input_list = [];

x_guess0 = 30;
num_trials = 100;
solver_flag = 2;

x_root = 37.879;
guess_list1 = linspace(x_root - 0.4, x_root + 0.4, num_trials);
guess_list2 = 0;
filter_list = [1e-9, 1e-2, 1e-9, 1e-2, 14];

% The quadratic's function value is the square of its root error
[p, k] = convergence_analysis(solver_flag, @fun, x_guess0, guess_list1, guess_list2, filter_list, x_root, 1e-26);

function [f_val,dfdx] = fun(x)
    global input_list;
    if isempty(input_list)
        input_list = x;
    else
        input_list(:,end+1) = x;
    end
    f_val = (x-37.879).^2;
    dfdx = 2*(x-37.879);
end
