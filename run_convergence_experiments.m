% runs convergence_analysis for all 4 solvers on 2 different functions
% and compiles predicted vs measured p and k into a comparison table
% Ensure all repo subfolders are in MATLAB's path regardless of working directory
root_dir = fileparts(mfilename('fullpath'));
addpath(genpath(root_dir));

funcs = {@test1, @test2};
names = {'test_func01', 'test_func02'};
solver_names = {'Bisection', 'Newton', 'Secant', 'fzero'};
x_guess0s = [0.5, 0.5];
num_trials = 100;

% Storage containers for table output
table_func = {};
table_solver = {};
table_p_pred = [];
table_p_meas = [];
table_k_pred = [];
table_k_meas = [];

for fi = 1:2
    fun = funcs{fi};
    x_guess0 = x_guess0s(fi);
    x_root = fzero(fun, x_guess0);
    fprintf('\n--- %s (root ~ %.6f) ---\n', names{fi}, x_root);
    
    for solver_flag = 1:4
        if solver_flag == 1
            % Bisection
            widths = linspace(0.05, 0.5, num_trials);
            ratios = linspace(0.3, 0.9, num_trials);
            guess_list1 = x_root - widths;
            guess_list2 = x_root + ratios .* widths;
            filter_list = [1e-15, 1e-2, 1e-14, 1e-2, 2];
            
            p_pred = 1;
            k_pred = 0.5;
            
        elseif solver_flag == 2
            % Newton
            guess_list1 = linspace(x_root - 0.4, x_root + 0.4, num_trials);
            guess_list2 = 0;
            filter_list = [1e-14, 1e-1, 1e-14, 1e-1, 1];
            
            [df, d2f] = approximate_derivative(fun, x_root);
            p_pred = 2;
            k_pred = abs(d2f / (2 * df));
            
        elseif solver_flag == 3
            % Secant
            widths = linspace(0.05, 0.5, num_trials);
            guess_list1 = x_root - widths;
            guess_list2 = x_root + 0.25 * widths;
            filter_list = [1e-14, 1e-2, 1e-14, 1e-2, 2];
            
            p_pred = 1.618;
            k_pred = NaN; % Rubric specifies no predicted k for secant
            
        elseif solver_flag == 4
            % fzero (Brent's method)
            guess_list1 = linspace(x_root - 0.4, x_root + 0.4, num_trials);
            guess_list2 = 0;
            filter_list = [1e-15, 1e-2, 1e-14, 1e-2, 2];
            
            p_pred = NaN; % Hybrid algorithm; no fixed theoretical rate
            k_pred = NaN;
        end
        
        [p, k] = convergence_analysis(solver_flag, fun, x_guess0, guess_list1, guess_list2, filter_list);
        
        % Append to table vectors
        table_func{end+1, 1} = names{fi};
        table_solver{end+1, 1} = solver_names{solver_flag};
        table_p_pred(end+1, 1) = p_pred;
        table_p_meas(end+1, 1) = p;
        table_k_pred(end+1, 1) = k_pred;
        table_k_meas(end+1, 1) = k;
        
        fprintf('solver %d: measured p=%.3f k=%.4f | predicted p=%.3f k=%.4f\n', ...
                solver_flag, p, k, p_pred, k_pred);
    end
end

% Construct and display comparison table
comparison_table = table(table_func, table_solver, table_p_pred, table_p_meas, table_k_pred, table_k_meas, ...
    'VariableNames', {'Function', 'Solver', 'Predicted_p', 'Measured_p', 'Predicted_k', 'Measured_k'});

fprintf('\n================== COMPARISON TABLE ==================\n');
disp(comparison_table);

% Run the Part 4 quadratic and sigmoid experiments
part4_results = run_part4_experiments();

function [fval,dfdx] = test1(x)
    fval = (x.^3)/100 - (x.^2)/8 + 2*x + 6*sin(x/2+6) -.7 - exp(x/6);
    dfdx = 3*(x.^2)/100 - 2*x/8 + 2 +(6/2)*cos(x/2+6) - exp(x/6)/6;
end

function [fval,dfdx] = test2(x)
    fval = exp(x) - 3*x;
    dfdx = exp(x) - 3;
end