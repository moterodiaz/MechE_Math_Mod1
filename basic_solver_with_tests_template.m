%template for testing your basic root finding implementations
function basic_solver_with_tests_template()
    xvals = linspace(-50,50,201);
    [yvals,~] = test_func01(xvals);

    hold on
    axis([-15,40,-50,80]);
    plot(xvals,yvals,'r','linewidth',2);
    plot(xvals,0*xvals,'k--','linewidth',1);
    xlabel('x'); ylabel('y'); title('Test Function 1');

    % %Newton's method example test
    % x0_guess = 2;
    % plot(x0_guess,test_func01(x0_guess),'bo','markerfacecolor','b','markersize',5);
    % 
    % x_sol = newton_solver(@test_func01,x0_guess);
    % plot(x_sol,test_func01(x_sol),'go','markerfacecolor','g','markersize',5);
    

    %Secant method example test
    % x0_guess = -5;
    % x1_guess = 2;
    % plot(x0_guess,test_func01(x0_guess),'bo','markerfacecolor','b','markersize',5);
    % plot(x1_guess,test_func01(x1_guess),'ko','markerfacecolor','k','markersize',5);
    % 
    % x_sol = secant_solver(@test_func01,x0_guess,x1_guess);
    % plot(x_sol,test_func01(x_sol),'go','markerfacecolor','g','markersize',5);

    
    %Bisection method example test
    x_left = -5;
    x_right = 2;
    dxtol = 0.0001;
    ftol = 0.0001;
    max_iter = 100;
    plot(x_left,test_func01(x_left),'bo','markerfacecolor','b','markersize',5);
    plot(x_right,test_func01(x_right),'ko','markerfacecolor','k','markersize',5);

    x_sol = bisection_solver(@test_func01,x_left,x_right,dxtol,ftol,max_iter);
    plot(x_sol,test_func01(x_sol),'go','markerfacecolor','g','markersize',5);

    fprintf('Approximate root: (x, y) = (%.6f, %.6f)\n', x_sol, test_func01(x_sol));
end


%Definition of the test function and its derivative (as a single function):
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
function x = newton_solver(fun,x0)
    x = x0+1; %this is just dummy code. replace this with your code
end

function x = secant_solver(fun,x0, x1)
    x = x0+1; %this is just dummy code. replace this with your code
end




