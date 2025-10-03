function Tasks = benchmark_TEST_MTSO(index)
%BENCHMARK function
%   Input
%   - index: the index number of problem set
%
%   Output:
%   - Tasks: benchmark problem set

%------------------------------- Reference --------------------------------
% @Article{Da2017CEC2017-MTSO,
%   author     = {Da, Bingshui and Ong, Yew-Soon and Feng, Liang and Qin, A Kai and Gupta, Abhishek and Zhu, Zexuan and Ting, Chuan-Kang and Tang, Ke and Yao, Xin},
%   journal    = {arXiv preprint arXiv:1706.03470},
%   title      = {Evolutionary Multitasking for Single-objective Continuous Optimization: Benchmark Problems, Performance Metric, and Baseline Results},
%   year       = {2017},
% }
%--------------------------------------------------------------------------

%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 Yanchi Li. You are free to use the MTO-Platform for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "MTO-Platform" and cite
% or footnote "https://github.com/intLyc/MTO-Platform"
%--------------------------------------------------------------------------

file_dir = './Problems/Multi-task/Rondom task pairs-MTSO/Tasks/';

switch (index)
    case 1 % complete intersection with high similarity, Griewank and Rastrigin
        load([file_dir, 'CI_H.mat']) % loading data from folder .\Tasks
        dim = 50;
        Tasks(1).Dim = dim; % dimensionality of Task 1
        Tasks(1).Fnc = @(x)Griewank(x, Rotation_Task1, GO_Task1, 0);
        Tasks(1).Lb = -100 * ones(1, dim); % Upper bound of Task 1
        Tasks(1).Ub = 100 * ones(1, dim); % Lower bound of Task 1

        Tasks(2).Dim = dim; % dimensionality of Task 2
        Tasks(2).Fnc = @(x)Rastrigin(x, Rotation_Task2, GO_Task2, 0);
        Tasks(2).Lb = -50 * ones(1, dim); % Upper bound of Task 2
        Tasks(2).Ub = 50 * ones(1, dim); % Lower bound of Task 2

        Tasks(3).Dim = dim; % dimensionality of Task 3
        Tasks(3).Fnc = @(x)Ackley(x, Rotation_Task1, GO_Task1, 0);
        Tasks(3).Lb = -50 * ones(1, dim); % Upper bound of Task 3
        Tasks(3).Ub = 50 * ones(1, dim); % Lower bound of Task 3

        Tasks(4).Dim = dim; % dimensionality of Task 4
        Tasks(4).Fnc = @(x)Rosenbrock(x, Rotation_Task2, GO_Task2, 0);
        Tasks(4).Lb = -50 * ones(1, dim); % Upper bound of Task 4
        Tasks(4).Ub = 50 * ones(1, dim); % Lower bound of Task 4

    case 2 % complete intersection with medium similarity, Ackley and Rastrigin
        load([file_dir, 'CI_H.mat'])
        dim = 50;
        
        Tasks(1).Dim = dim;
        Tasks(1).Fnc = @(x)Rastrigin(x, Rotation_Task1, GO_Task1, 0);
        Tasks(1).Lb = -50 * ones(1, dim);
        Tasks(1).Ub = 50 * ones(1, dim);

        Tasks(2).Dim = dim;
        Tasks(2).Fnc = @(x)Ackley(x, Rotation_Task2, GO_Task2, 0);
        Tasks(2).Lb = -50 * ones(1, dim);
        Tasks(2).Ub = 50 * ones(1, dim);

        Tasks(3).Dim = dim; % dimensionality of Task 3
        Tasks(3).Fnc = @(x)Griewank(x, Rotation_Task1, GO_Task1, 0);
        Tasks(3).Lb = -50 * ones(1, dim); % Upper bound of Task 3
        Tasks(3).Ub = 50 * ones(1, dim); % Lower bound of Task 3

        Tasks(4).Dim = dim; % dimensionality of Task 4
        Tasks(4).Fnc = @(x)Weierstrass(x, Rotation_Task2, GO_Task2, 0);
        Tasks(4).Lb = -50 * ones(1, dim); % Upper bound of Task 4
        Tasks(4).Ub = 50 * ones(1, dim); % Lower bound of Task 4

        Tasks(5).Dim = dim; % dimensionality of Task 5
        Tasks(5).Fnc = @(x)Ackley(x, Rotation_Task1, GO_Task1, 0);
        Tasks(5).Lb = -50 * ones(1, dim); % Upper bound of Task 5
        Tasks(5).Ub = 50 * ones(1, dim); % Lower bound of Task 5

        Tasks(6).Dim = dim; % dimensionality of Task 6
        Tasks(6).Fnc = @(x)Schwefel(x, Rotation_Task2, GO_Task2, 0);
        Tasks(6).Lb = -50 * ones(1, dim); % Upper bound of Task 6
        Tasks(6).Ub = 50 * ones(1, dim); % Lower bound of Task 6

    case 3 % complete intersection with low similarity, Ackley and Schwefel
        load([file_dir, 'CI_H.mat'])
        dim = 50;
        Tasks(1).Dim = dim;
        Tasks(1).Fnc = @(x)Schwefel(x, Rotation_Task1, GO_Task1, 0);
        Tasks(1).Lb = -50 * ones(1, dim);
        Tasks(1).Ub = 50 * ones(1, dim);

        Tasks(2).Dim = dim;
        Tasks(2).Fnc = @(x)Ackley(x, Rotation_Task2, GO_Task2, 0);
        Tasks(2).Lb = -500 * ones(1, dim);
        Tasks(2).Ub = 500 * ones(1, dim);

        Tasks(3).Dim = dim; % dimensionality of Task 3
        Tasks(3).Fnc = @(x)Griewank(x, Rotation_Task1, GO_Task1, 0);
        Tasks(3).Lb = -50 * ones(1, dim); % Upper bound of Task 3
        Tasks(3).Ub = 50 * ones(1, dim); % Lower bound of Task 3

        Tasks(4).Dim = dim; % dimensionality of Task 4
        Tasks(4).Fnc = @(x)Rosenbrock(x, Rotation_Task2, GO_Task2, 0);
        Tasks(4).Lb = -50 * ones(1, dim); % Upper bound of Task 4
        Tasks(4).Ub = 50 * ones(1, dim); % Lower bound of Task 4

        Tasks(5).Dim = dim; % dimensionality of Task 5
        Tasks(5).Fnc = @(x)Rastrigin(x, Rotation_Task1, GO_Task1, 0);
        Tasks(5).Lb = -50 * ones(1, dim); % Upper bound of Task 5
        Tasks(5).Ub = 50 * ones(1, dim); % Lower bound of Task 5

        Tasks(6).Dim = dim; % dimensionality of Task 6
        Tasks(6).Fnc = @(x)Weierstrass(x, Rotation_Task2, GO_Task2, 0);
        Tasks(6).Lb = -50 * ones(1, dim); % Upper bound of Task 6
        Tasks(6).Ub = 50 * ones(1, dim); % Lower bound of Task 6

        Tasks(7).Dim = dim; % dimensionality of Task 7
        Tasks(7).Fnc = @(x)Ackley(x, Rotation_Task1, GO_Task2, 0);
        Tasks(7).Lb = -50 * ones(1, dim); % Upper bound of Task 7
        Tasks(7).Ub = 50 * ones(1, dim); % Lower bound of Task 7

        Tasks(8).Dim = dim; % dimensionality of Task 8
        Tasks(8).Fnc = @(x)Rosenbrock(x, Rotation_Task2, GO_Task2, 0);
        Tasks(8).Lb = -50 * ones(1, dim); % Upper bound of Task 8
        Tasks(8).Ub = 50 * ones(1, dim); % Lower bound of Task 8

   
end
end
