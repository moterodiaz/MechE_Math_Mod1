function convergence_guesses(solver_flag, fun, x0_ref, function_name)
%solver_flag: 1 bisection 3 secant 4 fzero
%runs the solver for different initial guesses and plots success or failure
%Newton uses the existing newton_initial_guess function
% Ensure repository root and all subfolders are in MATLAB's path
repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repo_root));

    if ~ismember(solver_flag,[1,3,4])
        error('Use newton_initial_guess for the Newton experiment');
    end

    % ---------------------------------------------------------
    % Find the reference root
    % ---------------------------------------------------------

    target_root = fzero(fun,x0_ref);


    % ---------------------------------------------------------
    % Solver parameters
    % ---------------------------------------------------------

    dxtol = 0.0001;
    ftol = 0.0001;
    max_iter = 100;
    dxmax = 1000;
    root_tol = 0.001;
    % Tight fzero tolerance so its result passes the same residual check
    fzero_options = optimset('Display','off','TolX',1e-12);


    % ---------------------------------------------------------
    % Initial guesses
    % ---------------------------------------------------------

    num_trials = 201;
    x0_list = linspace(0,50,num_trials);
    x1_list = x0_list;

    if solver_flag == 1 || solver_flag == 3
        success_list = false(num_trials,num_trials);
    else
        success_list = false(1,num_trials);
    end

    root_list = NaN(size(success_list));
    exit_list = zeros(size(success_list));


    % ---------------------------------------------------------
    % Run the solver for every initial guess
    % ---------------------------------------------------------

    if solver_flag == 1 || solver_flag == 3

        for n = 1:num_trials
            for m = 1:num_trials

                x0 = x0_list(n);
                x1 = x1_list(m);

                if solver_flag == 1
                    [x_root,exit_flag] = bisection_solver( ...
                        fun,x0,x1,dxtol,ftol,max_iter);
                else
                    [x_root,exit_flag] = secant_solver( ...
                        fun,x0,x1,dxtol,ftol,max_iter,dxmax);
                end

                % Column is the horizontal guess, row is the vertical guess
                root_list(m,n) = x_root;
                exit_list(m,n) = exit_flag;
                success_list(m,n) = exit_flag == 1 && isfinite(x_root) && ...
                    abs(x_root-target_root) <= root_tol && abs(fun(x_root)) <= ftol;

            end
        end

    else

        for n = 1:num_trials

            x0 = x0_list(n);

            % One failed fzero guess should not stop the experiment
            try
                [x_root,~,exit_flag] = fzero(fun,x0,fzero_options);
            catch
                x_root = NaN;
                exit_flag = -1;
            end

            root_list(n) = x_root;
            exit_list(n) = exit_flag;
            success_list(n) = exit_flag == 1 && isfinite(x_root) && ...
                abs(x_root-target_root) <= root_tol && abs(fun(x_root)) <= ftol;

        end

    end


    % ---------------------------------------------------------
    % Plot convergence
    % ---------------------------------------------------------

    solver_names = {'Bisection','Newton','Secant','fzero'};
    name = solver_names{solver_flag};
    figure('Color','w','Position',[100,100,1100,850]);

    if solver_flag == 1 || solver_flag == 3

        imagesc(x0_list,x1_list,double(success_list),[0,1]);
        set(gca,'YDir','normal');
        colormap([1,0,0;0,0,1]);
        hold on;

        plot(NaN,NaN,'bs','markerfacecolor','b','DisplayName','Success');
        plot(NaN,NaN,'rs','markerfacecolor','r','DisplayName','Failure');
        xline(target_root,'k--','linewidth',1.5,'DisplayName','Root location');
        yline(target_root,'k--','linewidth',1.5,'HandleVisibility','off');
        plot(target_root,target_root,'go','markerfacecolor','g', ...
            'markersize',10,'DisplayName','(x_{root}, x_{root})');

        if solver_flag == 1
            xlabel('Initial value of x_L');
            ylabel('Initial value of x_R');
        else
            xlabel('Initial value of x_0');
            ylabel('Initial value of x_1');
        end

        axis square;

    else

        yvals = fun(x0_list);
        y_success = yvals;
        y_failure = yvals;
        y_success(~success_list) = NaN;
        y_failure(success_list) = NaN;

        plot(x0_list,y_success,'b.-','linewidth',2,'markersize',6, ...
            'DisplayName','Success');
        hold on;
        plot(x0_list,y_failure,'r.-','linewidth',2,'markersize',6, ...
            'DisplayName','Failure');
        yline(0,'k--','linewidth',1.5,'DisplayName','y = 0');
        plot(target_root,0,'go','markerfacecolor','g', ...
            'markersize',10,'DisplayName','Root');

        xlabel('x (initial guess x_0)');
        ylabel('f(x)');
        xlim([x0_list(1),x0_list(end)]);
        grid on;

    end

    title(sprintf('%s %s Guess Convergence',name,function_name), ...
        'FontSize',18,'Interpreter','none');
    set(gca,'FontSize',15);
    legend('Location','southoutside','Orientation','horizontal','FontSize',13);
    box on;
    hold off;

    fprintf('%s: %d of %d initial guesses converged to %.12f\n', ...
        name,nnz(success_list),numel(success_list),target_root);

end
