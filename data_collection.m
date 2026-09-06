function [e_n, e_np1, index_list] = data_collection(solver_flag, fun, x_root, guess_list1, guess_list2)
%solver_flag: 1 bisection, 2 newton, 3 secant, 4 fzero

dxtol = 1e-13;
ftol = 1e-13;
max_iter = 200;
dxmax = 1e10;

e_n = [];
e_np1 = [];
index_list = [];

my_recorder = input_recorder();
f_record = my_recorder.generate_recorder_fun(fun);

num_trials = length(guess_list1);

    for n = 1:num_trials

    if solver_flag == 1
        %bisection
        [~,~,iterates] = bisection_solver(fun, guess_list1(n), guess_list2(n), dxtol, ftol, max_iter);
    elseif solver_flag == 2
        %newton, recorder gets called twice per step so skip every other one
        my_recorder.clear_input_list();
        newton_solver(f_record, guess_list1(n), dxtol, ftol, max_iter, dxmax);
        input_list = my_recorder.get_input_list();
        iterates = input_list(1:2:end);
    elseif solver_flag == 3
        %secant
        my_recorder.clear_input_list();
        secant_solver(f_record, guess_list1(n), guess_list2(n), dxtol, ftol, max_iter, dxmax);
        iterates = my_recorder.get_input_list();
    elseif solver_flag == 4
        %fzero
        my_recorder.clear_input_list();
        fzero(f_record, guess_list1(n));
        iterates = my_recorder.get_input_list();
    end

    %need at least 2 points to make a pair
    if length(iterates) >= 2
        e_n = [e_n, abs(iterates(1:end-1) - x_root)];
        e_np1 = [e_np1, abs(iterates(2:end) - x_root)];
        index_list = [index_list, 1:length(iterates)-1];
    end

    end

end
