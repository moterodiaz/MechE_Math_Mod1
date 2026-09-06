function [p, k] = convergence_analysis(solver_flag, fun, x_guess0, guess_list1, guess_list2, filter_list)
%solver_flag: 1 bisection 2 newton 3 secant 4 fzero
%filter_list = [e_n min, e_n max, e_np1 min, e_np1 max, min iteration]

x_root = fzero(fun, x_guess0);

[e_n, e_np1, index_list] = data_collection(solver_flag, fun, x_root, guess_list1, guess_list2);

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

figure;
loglog(e_n, e_np1, 'ro','markerfacecolor','r','markersize',1);
hold on
loglog(x_regression, y_regression, 'bo','markerfacecolor','b','markersize',2);

fit_line_x = 10.^[-16:.01:1];
fit_line_y = k*fit_line_x.^p;
loglog(fit_line_x, fit_line_y, 'k-','linewidth',2);

xlabel('e_n'); ylabel('e_{n+1}');
title(name);
grid on

fprintf('%s: root=%.10f p=%.4f k=%.4f\n', name, x_root, p, k);

end
