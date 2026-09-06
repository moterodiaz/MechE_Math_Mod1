%finite difference approximation for f'(x) and f''(x)
%INPUTS:
%   fun: the function we want to differentiate (fun(x) should output just f)
%   x: the point to differentiate at
%OUTPUTS:
%   dfdx: approximation of fun'(x)
%   d2fdx2: approximation of fun''(x)
function [dfdx,d2fdx2] = approximate_derivative(fun,x)
    delta_x = 1e-6;

    f_left = fun(x-delta_x);
    f_0 = fun(x);
    f_right = fun(x+delta_x);

    dfdx = (f_right-f_left)/(2*delta_x);
    d2fdx2 = (f_right-2*f_0+f_left)/(delta_x^2);
end
