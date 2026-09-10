function [x_s, x_f, f_s, f_f] = newton_initial_guess()
% Track which initial guesses reach the root
% and the function values at each initial guess.

% Prompt user for function inputs
min_guess = input('Enter minimum initial guess: ');
max_guess = input('Enter maximum initial guess: ');
trials = input('Enter number of trials: ');

% Create an array of initial guesses
guess_range = linspace(min_guess, max_guess, trials);

dxtol = 0.0001;     % Distance tolerance
ftol = 0.0001;      % Function tolerance
max_iter = 100;     % Maximum number of iterations
dxmax = 1000;       % Maximum allowed Newton step

% Declare function used for testing
function [f_val,dfdx] = test_function03(x)
    a = 27.3;
    b = 2;
    c = 8.3;
    d = -3;

    H = exp((x-a)/b);
    dH = H/b;
    L = 1+H;
    dL = dH;

    f_val = c*H./L+d;
    dfdx = c*(L.*dH-H.*dL)./(L.^2);
end

% Create arrays for successful and failed initial guesses
x_s = NaN(1,trials);
x_f = NaN(1,trials);

% Create arrays for function values
f_s = NaN(1,trials);
f_f = NaN(1,trials);

% Test each initial guess
for i = 1:trials

    % Current initial guess
    x0 = guess_range(i);

    % Function value at the initial guess
    [f_value, ~] = test_function03(x0);

    % Run Newton's method
    [x, exit_flag] = newton_solver(@test_function03, ...
        x0, dxtol, ftol, max_iter, dxmax);

    % Separate successful and failed guesses
    if exit_flag == 1
        x_s(i) = x0;
        f_s(i) = f_value;
    else
        x_f(i) = x0;
        f_f(i) = f_value;
    end

end

% Plot successful initial guesses
plot(x_s, f_s, 'o', 'Color', 'Green');
xlabel('Trial');
ylabel('Successful Initial Guess');
grid on;
xlim([min_guess max_guess])
hold on;
plot(x_f, f_f, 'o', 'Color', 'red');

end