function part4_results = run_part4_experiments()
% runs the Part 4 quadratic tests and initial guess experiments
% and compiles predicted vs measured p and k into a comparison table
root_dir = fileparts(mfilename('fullpath'));
addpath(genpath(root_dir));

solver_names = {'Bisection','Newton','Secant','fzero'};
num_trials = 100;

% Storage containers for table output
table_solver = {};
table_p_pred = [];
table_p_meas = [];
table_k_pred = [];
table_k_meas = [];

for solver_flag = 2:3

    if solver_flag == 2
        % Use the existing Newton test for the quadratic at its minimum
        testnewton;
        p_pred = 1;
        k_pred = 0.5;
    else
        % Secant, use different guesses on the same side of the root
        x_root = 37.879;
        x_guess0 = 30;
        widths = linspace(0.05,0.5,num_trials);
        guess_list1 = x_root + widths;
        guess_list2 = x_root + 0.6*widths;
        filter_list = [1e-9,1e-2,1e-9,1e-2,14];

        [p,k] = convergence_analysis(solver_flag,@test_function02,x_guess0, ...
            guess_list1,guess_list2,filter_list,x_root,1e-26);
        p_pred = 1;
        k_pred = NaN; % No predicted k is required for secant
    end

    % Append to table vectors
    table_solver{end+1,1} = solver_names{solver_flag};
    table_p_pred(end+1,1) = p_pred;
    table_p_meas(end+1,1) = p;
    table_k_pred(end+1,1) = k_pred;
    table_k_meas(end+1,1) = k;

end

% Construct and display comparison table
comparison_table = table(table_solver,table_p_pred,table_p_meas,table_k_pred,table_k_meas, ...
    'VariableNames',{'Solver','Predicted_p','Measured_p','Predicted_k','Measured_k'});
comparison_table.Predicted_k = string(comparison_table.Predicted_k);
comparison_table.Predicted_k(strcmp(table_solver,'Secant')) = "Not required";

fprintf('\n================== PART 4 COMPARISON TABLE ==================\n');
disp(comparison_table);

% Use the existing initial guess experiment for Newton
newton_initial_guess(0,50,201);

% Run the other initial guess experiments on the same sigmoid
for solver_flag = [1,3,4]
    convergence_guesses(solver_flag,@test_function03,27.3,'Sigmoid Function');
end
part4_results.comparison_table = comparison_table;

end

%Quadratic function with root at the minimum
function [f_val,dfdx] = test_function02(x)
    f_val = (x-37.879).^2;
    dfdx = 2*(x-37.879);
end

%Example sigmoid function
function [f_val,dfdx] = test_function03(x)
    a = 27.3; b = 2; c = 8.3; d = -3;
    H = exp((x-a)/b);
    dH = H/b;
    L = 1+H;
    dL = dH;
    f_val = c*H./L+d;
    dfdx = c*(L.*dH-H.*dL)./(L.^2);
end
