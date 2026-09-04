function convergence_secant()

    % ---------------------------------------------------------
    % Find the reference root
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

    dxtol = 1e-14;
    ftol = 1e-14;
    max_iter = 100;
    dxmax = 10;


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
    % Run many Secant experiments
    % ---------------------------------------------------------

    for n = 1:num_trials

        % Choose different starting guesses

        width = 0.05 + ...
                0.45*(n-1)/(num_trials-1);

        x0 = target_root - width;
        x1 = target_root + 0.25*width;


        % Clear recorder

        my_recorder.clear_input_list();


        % Run Secant method

        [x_root,exit_flag] = secant_solver( ...
            f_record,...
            x0,...
            x1,...
            dxtol,...
            ftol,...
            max_iter,...
            dxmax);


        % Get recorded guesses

        input_list = my_recorder.get_input_list();


        % -----------------------------------------------------
        % Create convergence pairs
        % -----------------------------------------------------

        if length(input_list) >= 3

            current = input_list(1:end-1);
            next = input_list(2:end);


            % Calculate errors for this trial

            current_error = ...
                abs(current-target_root);

            next_error = ...
                abs(next-target_root);


            % Keep valid points

            valid = ...
                isfinite(current_error) & ...
                isfinite(next_error) & ...
                current_error > 0 & ...
                next_error > 0;


            current_error = current_error(valid);
            next_error = next_error(valid);


            % Add this trial's data

            x_current_list = ...
                [x_current_list current_error];

            x_next_list = ...
                [x_next_list next_error];

        end

    end


    % ---------------------------------------------------------
    % Remove extreme numerical errors
    % ---------------------------------------------------------

    valid = ...
        x_current_list > 1e-14 & ...
        x_next_list > 1e-14 & ...
        x_current_list < 10 & ...
        x_next_list < 10;

    x_current_list = ...
        x_current_list(valid);

    x_next_list = ...
        x_next_list(valid);


    % ---------------------------------------------------------
    % Display results
    % ---------------------------------------------------------

    fprintf('Target root = %.15f\n',target_root);

    fprintf('Number of convergence points = %d\n',...
        length(x_current_list));


    % ---------------------------------------------------------
    % Estimate convergence order
    % ---------------------------------------------------------

    if length(x_current_list) >= 3

        log_x = log10(x_current_list);
        log_y = log10(x_next_list);

        p = polyfit(log_x,log_y,1);

        fprintf('Estimated convergence order = %.4f\n',p(1));

    end


    % ---------------------------------------------------------
    % Plot convergence
    % ---------------------------------------------------------

    figure;

    loglog(x_current_list,...
           x_next_list,...
           'ro',...
           'markersize',3);

    xlabel('|e_n|');
    ylabel('|e_{n+1}|');

    title('Secant Method Error Convergence');

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