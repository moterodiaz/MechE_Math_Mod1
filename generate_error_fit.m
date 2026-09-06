%fits a line to log(error) data to get p (order) and k (constant)
%from the convergence relation e_{n+1} = k*(e_n)^p
%INPUTS:
%   x_regression: e_n values
%   y_regression: e_{n+1} values
%OUTPUTS:
%   p: order of convergence (slope of the fit)
%   k: convergence constant (intercept of the fit, undone from log)
function [p,k] = generate_error_fit(x_regression,y_regression)
    Y = log(y_regression)';
    X1 = log(x_regression)';
    X2 = ones(length(X1),1);

    coeff_vec = regress(Y,[X1,X2]);

    p = coeff_vec(1);
    k = exp(coeff_vec(2));
end
