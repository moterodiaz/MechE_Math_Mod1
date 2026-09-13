%fits a line to log(error) data to get p (order) and k (constant)
%from the convergence relation e_{n+1} = k*(e_n)^p
%INPUTS:
%   x_regression: e_n values
%   y_regression: e_{n+1} values
%OUTPUTS:
%   p: order of convergence (slope of the fit)
%   k: convergence constant (intercept of the fit, undone from log)
function [p,k] = generate_error_fit(x_regression,y_regression)
    x_regression = x_regression(:);
    y_regression = y_regression(:);
    valid = isfinite(x_regression) & isfinite(y_regression) & ...
            x_regression > 0 & y_regression > 0;
    x_regression = x_regression(valid);
    y_regression = y_regression(valid);

    if length(x_regression) < 3 || length(unique(x_regression)) < 2
        error('Not enough error pairs to fit p and k');
    end

    Y = log(y_regression);
    X1 = log(x_regression);
    X2 = ones(length(X1),1);

    coeff_vec = [X1,X2]\Y;

    p = coeff_vec(1);
    k = exp(coeff_vec(2));
end
