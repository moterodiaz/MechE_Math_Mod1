%runs convergence_analysis for all 4 solvers on 2 different functions
%and prints predicted vs measured p and k

funcs = {@test1, @test2};
names = {'test_func01','test_func02'};
x_guess0s = [0.5, 0.5];

filter_list = [1e-15, 1e-2, 1e-14, 1e-2, 2];
num_trials = 100;

for fi = 1:2

    fun = funcs{fi};
    x_guess0 = x_guess0s(fi);
    x_root = fzero(fun, x_guess0); %just used to build guess ranges around

    fprintf('\n--- %s (root ~ %.6f) ---\n', names{fi}, x_root);

    for solver_flag = 1:4

        if solver_flag == 1
            widths = linspace(0.05,0.5,num_trials);
            ratios = linspace(0.3,0.9,num_trials); %varying, not a fixed ratio
            guess_list1 = x_root - widths;
            guess_list2 = x_root + ratios.*widths;
        elseif solver_flag == 2
            guess_list1 = linspace(x_root-0.4, x_root+0.4, num_trials);
            guess_list2 = 0;
        elseif solver_flag == 3
            widths = linspace(0.05,0.5,num_trials);
            guess_list1 = x_root - widths;
            guess_list2 = x_root + 0.25*widths;
        elseif solver_flag == 4
            guess_list1 = linspace(x_root-0.4, x_root+0.4, num_trials);
            guess_list2 = 0;
        end

        [p,k] = convergence_analysis(solver_flag, fun, x_guess0, guess_list1, guess_list2, filter_list);

        if solver_flag == 1
            p_pred = 1;
            k_pred = 0.5;
        elseif solver_flag == 2
            [df,d2f] = approximate_derivative(fun, x_root);
            p_pred = 2;
            k_pred = abs(d2f/(2*df));
        elseif solver_flag == 3
            p_pred = 1.618;
            k_pred = NaN; %no simple formula for secant
        else
            p_pred = NaN; %fzero isn't just one method, no single prediction
            k_pred = NaN;
        end

        fprintf('solver %d: measured p=%.3f k=%.4f | predicted p=%.3f k=%.4f\n', solver_flag, p, k, p_pred, k_pred);

    end
end


function [fval,dfdx] = test1(x)
    fval = (x.^3)/100 - (x.^2)/8 + 2*x + 6*sin(x/2+6) -.7 - exp(x/6);
    dfdx = 3*(x.^2)/100 - 2*x/8 + 2 +(6/2)*cos(x/2+6) - exp(x/6)/6;
end

function [fval,dfdx] = test2(x)
    fval = exp(x) - 3*x;
    dfdx = exp(x) - 3;
end
