function convergence_bisection()

    % ---------------------------------------------------------
    % Find the root
    % ---------------------------------------------------------

    x0_ref = 0.5;
    target_root = fzero(@test_func01,x0_ref);


    % ---------------------------------------------------------
    % Create input recorder
    % ---------------------------------------------------------

    my_recorder = input_recorder();

    f_record = ...
        my_recorder.generate_recorder_fun(@test_func01);


    % ---------------------------------------------------------
    % Solver parameters
    % ---------------------------------------------------------

    dxtol = 1e-12;
    ftol = 1e-12;
    max_iter = 200;


    % ---------------------------------------------------------
    % Number of trials
    % ---------------------------------------------------------

    num_trials = 100;


    % ---------------------------------------------------------
    % Lists for convergence data
    % ---------------------------------------------------------

    x_current_list = [];
    x_next_list = [];


    % ---------------------------------------------------------
    % Run many bisection experiments
    % ---------------------------------------------------------

    for n = 1:num_trials

        width = 0.05 + ...
                (0.45)*(n-1)/(num_trials-1);


        % Use an asymmetric bracket so that the
        % first midpoint is not exactly the root

        x_left = target_root - width;
        x_right = target_root + 0.6*width;


        % Check that the interval brackets the root

        if sign(test_func01(x_left)) == ...
           sign(test_func01(x_right))

            continue;

        end


        % Clear recorder

        my_recorder.clear_input_list();


        % Run bisection

        x_root = bisection_solver( ...
            f_record,...
            x_left,...
            x_right,...
            dxtol,...
            ftol,...
            max_iter);


        % Get recorded inputs

        input_list = my_recorder.get_input_list();


        % Extract midpoint values
        %
        % The solver records:
        % x_left, x_right, c1, x_left, c2, x_left, ...

        if length(input_list) >= 3

            midpoint_list = input_list(3:2:end);


            % Create convergence pairs

            if length(midpoint_list) >= 2

                x_current_list = ...
                    [x_current_list ...
                     midpoint_list(1:end-1)];

                x_next_list = ...
                    [x_next_list ...
                     midpoint_list(2:end)];

            end

        end

    end


    % ---------------------------------------------------------
    % Calculate errors
    % ---------------------------------------------------------

    abs_error_current = ...
        abs(x_current_list-target_root);

    abs_error_next = ...
        abs(x_next_list-target_root);


    % ---------------------------------------------------------
    % Remove invalid data
    % ---------------------------------------------------------

    valid = ...
        isfinite(abs_error_current) & ...
        isfinite(abs_error_next) & ...
        abs_error_current > 0 & ...
        abs_error_next > 0;

    abs_error_current = ...
        abs_error_current(valid);

    abs_error_next = ...
        abs_error_next(valid);


    % ---------------------------------------------------------
    % Remove significantly erroneous data points
    % ---------------------------------------------------------

    % Work in log space since the plot is log-log

    log_x = log10(abs_error_current);
    log_y = log10(abs_error_next);


    % Fit the main convergence trend

    p = polyfit(log_x,log_y,1);


    % Calculate distance from the fitted line

    log_y_fit = polyval(p,log_x);

    residual = abs(log_y-log_y_fit);


    % Keep points close to the main convergence trend

    threshold = 0.25;

    valid = residual < threshold;


    abs_error_current = ...
        abs_error_current(valid);

    abs_error_next = ...
        abs_error_next(valid);


    % ---------------------------------------------------------
    % Display results
    % ---------------------------------------------------------

    fprintf('Target root = %.15f\n',target_root);

    fprintf('Number of convergence points after filtering = %d\n',...
        length(abs_error_current));


    % ---------------------------------------------------------
    % Plot convergence
    % ---------------------------------------------------------

    figure;

    loglog(abs_error_current,...
           abs_error_next,...
           'ro',...
           'markerfacecolor','r',...
           'markersize',3);

    xlabel('\epsilon_n');
    ylabel('\epsilon_{n+1}');

    title('Bisection Method Error Convergence');

    grid on;

end


% =============================================================
% Test function
% =============================================================

function [fval,dfdx] = test_func01(x)

    fval = (x.^3)/100 ...
        - (x.^2)/8 ...
        + 2*x ...
        + 6*sin(x/2+6) ...
        -.7 ...
        - exp(x/6);

    dfdx = 3*(x.^2)/100 ...
        - 2*x/8 ...
        + 2 ...
        +(6/2)*cos(x/2+6) ...
        - exp(x/6)/6;

end