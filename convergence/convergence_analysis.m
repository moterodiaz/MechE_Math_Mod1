function [p, k] = convergence_analysis(solver_flag, fun, x_guess0, guess_list1, guess_list2, filter_list, x_root, ftol)
%solver_flag: 1 bisection 2 newton 3 secant 4 fzero
%filter_list = [e_n min, e_n max, e_np1 min, e_np1 max, min iteration]
% Ensure repository root and all subfolders are in MATLAB's path
repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repo_root));

if nargin < 7 || isempty(x_root) || x_root == 0
    x_root = fzero(fun, x_guess0);
end

if nargin < 8
    ftol = 1e-13;
end

[e_n, e_np1, index_list] = data_collection(solver_flag, fun, x_root, guess_list1, guess_list2, ftol);

%clean out the junk (too small = machine precision noise, too big = not converged yet)
x_regression = [];
y_regression = [];
for n = 1:length(index_list)
    if e_n(n) > filter_list(1) && e_n(n) < filter_list(2) && e_np1(n) > filter_list(3) && e_np1(n) < filter_list(4) && index_list(n) > filter_list(5)
        x_regression(end+1) = e_n(n);
        y_regression(end+1) = e_np1(n);
    end
end

[p,k] = generate_error_fit(x_regression, y_regression);

if solver_flag == 1
    name = 'Bisection';
elseif solver_flag == 2
    name = 'Newton';
elseif solver_flag == 3
    name = 'Secant';
elseif solver_flag == 4
    name = 'fzero';
end

% ---------------------------------------------------------
% Plotting and formatting
% ---------------------------------------------------------
figure;
loglog(e_n, e_np1, 'ro', 'markerfacecolor', 'r', 'markersize', 2, 'DisplayName', 'Unfiltered Data');
hold on;
loglog(x_regression, y_regression, 'bo', 'markerfacecolor', 'b', 'markersize', 3, 'DisplayName', 'Filtered Data');

% Bound fit line strictly to the filtered linear regime to prevent y-axis blowout
if ~isempty(x_regression)
    fit_line_x = logspace(log10(min(x_regression)), log10(max(x_regression)), 100);
    fit_line_y = k * (fit_line_x .^ p);
    loglog(fit_line_x, fit_line_y, 'k-', 'linewidth', 2, ...
           'DisplayName', sprintf('Fit Line (p = %.2f, k = %.2f)', p, k));
end
hold off;

% Labels, formatting, and legibility sizing
xlabel('e_n', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('e_{n+1}', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('%s Error Convergence Rate', name), 'FontSize', 13, 'FontWeight', 'bold');

set(gca, 'FontSize', 11);
xlim([1e-16, 1e1]);
ylim([1e-16, 1e1]);
grid on;
grid minor;

% Northwest corner is always unoccupied for convergent sequences (e_{n+1} < e_n)
legend('Location', 'northwest', 'FontSize', 10);

fprintf('%s: root=%.10f p=%.4f k=%.4f\n', name, x_root, p, k);
end
