% Convergence experiment for Newton's method
function convergence_newton()

    % ---------------------------------------------------------
    % Find the root we want to study
    % ---------------------------------------------------------

    x0_ref = 0.5;

    % Use MATLAB's fzero to find an accurate reference root
    target_root = fzero(@test_func01,x0_ref);


    % ---------------------------------------------------------
    % Create input recorder
    % ---------------------------------------------------------

    my_recorder = input_recorder();

    % Create a version of test_func01 that records every input
    f_record = my_recorder.generate_recorder_fun(@test_func01);


    % ---------------------------------------------------------
    % Number of trials
    % ---------------------------------------------------------

    num_trials = 1000;


    % ---------------------------------------------------------
    % Newton solver parameters
    % ---------------------------------------------------------

    dxtol = 1e-12;
    ftol = 1e-12;
    max_iter = 200;
    dxmax = 1e10;


    % ---------------------------------------------------------
    % Initial guesses
    % ---------------------------------------------------------

    x0_list = linspace(x0_ref-2,x0_ref+2,num_trials);


    % ---------------------------------------------------------
    % Lists to store convergence data
    % ---------------------------------------------------------

    x_current_list = [];
    x_next_list = [];
    index_list = [];


    % ---------------------------------------------------------
    % Run Newton's method for every initial guess
    % ---------------------------------------------------------

    for n = 1:num_trials

        % Get initial guess for this trial
        x0 = x0_list(n);

        % Clear recorder
        my_recorder.clear_input_list();

        % Run Newton's method
        [x_root,exit_flag] = newton_solver( ...
            f_record,x0,dxtol,ftol,max_iter,dxmax);

        % Get all recorded function inputs
        input_list = my_recorder.get_input_list();


        % -----------------------------------------------------
        % Extract the actual Newton iterates
        % -----------------------------------------------------
        %
        % With your newton_solver and input_recorder,
        % the recorded sequence is:
        %
        % x0, x1, x1, x2, x2, x3, x3, ...
        %
        % Therefore, every other value gives:
        %
        % x0, x1, x2, x3, ...
        %

        if length(input_list) >= 3

            iterates = input_list(1:2:end);


            % Need at least two iterates to calculate
            % an error pair
            if length(iterates) >= 2

                % Store x_n
                x_current_list = [x_current_list, ...
                                  iterates(1:end-1)];

                % Store x_(n+1)
                x_next_list = [x_next_list, ...
                               iterates(2:end)];

                % Store iteration number
                index_list = [index_list, ...
                              1:length(iterates)-1];

            end

        end

    end


    % ---------------------------------------------------------
    % Calculate errors
    % ---------------------------------------------------------

    abs_error_current = abs(x_current_list - target_root);
    abs_error_next = abs(x_next_list - target_root);


    % ---------------------------------------------------------
    % Remove invalid convergence data
    % ---------------------------------------------------------

    valid = abs_error_current > 0 & ...
            abs_error_next > 0 & ...
            isfinite(abs_error_current) & ...
            isfinite(abs_error_next);

    abs_error_current = abs_error_current(valid);
    abs_error_next = abs_error_next(valid);


    % ---------------------------------------------------------
    % Plot convergence
    % ---------------------------------------------------------

    figure;

    loglog(abs_error_current,abs_error_next, ...
        'ro','markerfacecolor','r','markersize',2);

    xlabel('|e_n|');
    ylabel('|e_{n+1}|');

    title('Newton''s Method Error Convergence');

    grid on;


    % ---------------------------------------------------------
    % Display results
    % ---------------------------------------------------------

    fprintf('Target root = %.12f\n',target_root);
    fprintf('Number of convergence data points = %d\n', ...
        length(abs_error_current));

end


% =============================================================
% Test function and derivative
% =============================================================

function [fval,dfdx] = test_func01(x)

    % Function
    fval = (x.^3)/100 ...
         - (x.^2)/8 ...
         + 2*x ...
         + 6*sin(x/2+6) ...
         - .7 ...
         - exp(x/6);

    % Derivative
    dfdx = 3*(x.^2)/100 ...
         - 2*x/8 ...
         + 2 ...
         + (6/2)*cos(x/2+6) ...
         - exp(x/6)/6;

end