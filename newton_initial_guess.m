function [x_s, x_f, f_s, f_f] = newton_initial_guess(min_guess, max_guess, trials)
% Track which initial guesses reach the root
% and the function values at each initial guess.

% Ensure all repo subfolders are in MATLAB's path
root_dir = fileparts(mfilename('fullpath'));
addpath(genpath(root_dir));

% Prompt user for function inputs
if nargin < 1
    min_guess = input('Enter minimum initial guess: ');
end
if nargin < 2
    max_guess = input('Enter maximum initial guess: ');
end
if nargin < 3
    trials = input('Enter number of trials: ');
end

% Create an array of initial guesses
guess_range = linspace(min_guess, max_guess, trials);

dxtol = 0.0001;     % Distance tolerance
ftol = 0.0001;      % Function tolerance
max_iter = 100;     % Maximum number of iterations
dxmax = 1000;       % Maximum allowed Newton step
target_root = 27.3 + 2*log(3/5.3);
root_tol = 0.001;

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
root_list = NaN(1,trials);
exit_list = zeros(1,trials);

% Test each initial guess
for i = 1:trials

    % Current initial guess
    x0 = guess_range(i);

    % Function value at the initial guess
    [f_value, ~] = test_function03(x0);

    % Run Newton's method
    [x, exit_flag] = newton_solver(@test_function03, ...
        x0, dxtol, ftol, max_iter, dxmax);
    root_list(i) = x;
    exit_list(i) = exit_flag;

    % Separate successful and failed guesses
    if exit_flag == 1 && isfinite(x) && ...
            abs(x-target_root) <= root_tol && abs(test_function03(x)) <= ftol
        x_s(i) = x0;
        f_s(i) = f_value;
    else
        x_f(i) = x0;
        f_f(i) = f_value;
    end

end

% Plot successful initial guesses
figure('Color','w','Position',[100,100,1100,750]);
plot(x_s, f_s, '.-', 'Color', 'blue', 'linewidth',2, ...
    'DisplayName','Success');
xlabel('x (initial guess x_0)');
ylabel('f(x)');
grid on;
xlim([min_guess max_guess])
hold on;
plot(x_f, f_f, '.-', 'Color', 'red', 'linewidth',2, ...
    'DisplayName','Failure');
yline(0,'k--','linewidth',1.5,'DisplayName','y = 0');
plot(target_root,0,'go','markerfacecolor','g','markersize',10, ...
    'DisplayName','Root');
title('Newton''s Method Sigmoid Function Guess Convergence','FontSize',18);
set(gca,'FontSize',15);
legend('Location','southoutside','Orientation','horizontal','FontSize',13);
box on;
hold off;

% Display summary
success_list = ~isnan(x_s);

fprintf('Newton: %d of %d initial guesses converged to %.12f\n', ...
    nnz(success_list),trials,target_root);

end
