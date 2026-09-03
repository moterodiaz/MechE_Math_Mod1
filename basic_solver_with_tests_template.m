%template for testing your basic root finding implementations
function basic_solver_with_tests_template()
    xvals = linspace(-50,50,201);
    [yvals,~] = test_func01(xvals);

    hold on
    axis([-15,40,-50,80]);
    plot(xvals,yvals,'r','linewidth',2);
    plot(xvals,0*xvals,'k--','linewidth',1);
    xlabel('x'); ylabel('y'); title('Test Function 1');

    %Newton's method example test
    x0_guess = 2;
    dxtol = 0.0001;
    ftol = 0.0001;
    max_iter = 100;
    dxmax = 1000;

    x_newton = newton_solver(@test_func01,x0_guess,dxtol,ftol,max_iter,dxmax);
    y_newton = test_func01(x_newton);
    plot(x_newton,y_newton,'go','markerfacecolor','g','markersize',5);
    

    %Secant method example test
    x0_guess = -5;
    x1_guess = 2;
    plot(x0_guess,test_func01(x0_guess),'bo','markerfacecolor','b','markersize',5);
    plot(x1_guess,test_func01(x1_guess),'ko','markerfacecolor','k','markersize',5);

    x_sol = secant_solver(@test_func01,x0_guess,x1_guess,dxtol,ftol,max_iter,dxmax);
    plot(x_sol,test_func01(x_sol),'go','markerfacecolor','g','markersize',5);

    
    %Bisection method example test
    x_left = -5;
    x_right = 2;
    dxtol = 0.0001;
    ftol = 0.0001;
    max_iter = 100;

    x_bisection = bisection_solver(@test_func01,x_left,x_right,dxtol,ftol,max_iter);
    y_bisection = test_func01(x_bisection);
    plot(x_bisection,y_bisection,'go','markerfacecolor','g','markersize',5);

    % Display both coordinate pairs
    fprintf('Newton Method:    (x, y) = (%.6f, %.6f)\n', x_newton, y_newton);
    fprintf('Bisection Method: (x, y) = (%.6f, %.6f)\n', x_bisection, y_bisection);
    fprintf('Secant Method: (x, y) = (%.6f, %.6f)\n', x_sol, test_func01(x_sol));

end


% Definition of the test function and its derivative (as a single function):
%This definition uses the function keyword
%when passing this function as an argument to a solver,
%you'll need to use the handle operator
%ex. solver(@test_func01,x_guess)
function [fval,dfdx] = test_func01(x)
    fval = (x.^3)/100 - (x.^2)/8 + 2*x + 6*sin(x/2+6) -.7 - exp(x/6);
    dfdx = 3*(x.^2)/100 - 2*x/8 + 2 +(6/2)*cos(x/2+6) - exp(x/6)/6;
end


function x = bisection_solver(fun,x_left,x_right,dxtol,ftol,max_iter) 
    % Initialize exit flag 
    exit_flag = 0; 
 
    if sign(fun(x_left)) == sign(fun(x_right)) 
        disp('Root Does Not Necessarily Exist Between Bounds') 
        return; 
    end 
 
    % Loop function max amount of times if necessary 
    for i = 1 : max_iter 
        
    % Calculate Midpoint 
        c = (x_left + x_right)/2; 
        x = c; 
 
    % Calculate function at midpoint 
    f_c = fun(c); 
     
    % Termination threshold --> Diff between estimate and zero 
    if abs(f_c) <= ftol 
         exit_flag = 1; 
         return; 
    end 
 
    % Create new bounds using midpoint 
         
    % If the left bound and midpoint have different sign values, midpoint 
    % becomes new right bound 
    % If the right bound and midpoint have different sign values, midpoint 
    % becomes new left bound 
    if sign(fun(x_left)) ~= sign(f_c) 
        x_right = c; 
    else 
        x_left = c; 
    end 
     
    % Termination threshold --> Diff between left and right bounds 
        if abs(x_left - x_right) <= dxtol 
            exit_flag = 1; 
            disp('X termination criteria met'); 
            return; 
        end 

    end
end

%Note that fun(x) should output [f,dfdx], where dfdx is the derivative of f 
function x = newton_solver(fun,x0,dxtol,ftol,max_iter,dxmax)

    % Initial root guess
    x_i = x0;

    for n = 1:max_iter

        % Get function and derivative
        [f_i,df_i] = fun(x_i);

        % Check for division by zero
        if df_i == 0
            x = x_i;
            return;
        end

        % Find the next guess
        x = x_i - f_i/df_i;

        % Check for excessively large change
        if abs(x - x_i) > dxmax
            return;
        end

        % Check function value at new guess
        f_x = fun(x);

        if abs(f_x) <= ftol
            return;
        end

        % Check how much the root guess changes
        if abs(x - x_i) <= dxtol
            return;
        end

        % Update guess
        x_i = x;
    end
    end
   
function x = secant_solver(fun, x0_guess, x1_guess, dxtol, ftol, max_iter, dxmax)

    x_first = x0_guess;
    x_second = x1_guess;
    denominator_tol = 1e-14;

    f_first = fun(x_first);
    f_second = fun(x_second);

    % check if guess is a root
    if abs(f_first) < ftol
        x = x_first;
        exit_flag = 1;
        return
    elseif abs(f_second) < ftol
        x = x_second;
        exit_flag = 1;
        return
    end

    %iterations until exit
    for iter = 1:max_iter
        denominator = f_second - f_first;

        % make sure denominator passes tolerance
        if ~isfinite(denominator) || abs(denominator) < denominator_tol
            x = x_second;
            exit_flag = -2;
            return
        end

        %secant method algorithm
        dx = -f_second * (x_second - x_first) / denominator;

        % terminate if dx is too large
        if ~isfinite(dx) || abs(dx) > dxmax
            x = x_second;
            exit_flag = -3;
            return
        end

        %secant update
        x_new = x_second + dx;
        f_new = fun(x_new);

        % early termination if near solution
        if abs(dx) < dxtol || abs(f_new) < ftol
            x = x_new;
            exit_flag = 1;
            return
        end

        % shift secant points forward one iteration on update
        x_first = x_second;
        f_first = f_second;

        x_second = x_new;
        f_second = f_new;
    end

    x = x_second;
    exit_flag = 0;
end




