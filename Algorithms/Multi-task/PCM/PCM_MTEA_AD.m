classdef PCM_MTEA_AD < Algorithm
% <Multi-task/Many-task> <Single-objective> <None/Constrained>

%------------------------------- Reference --------------------------------
% @Article{Wang2021MTEA-AD,
%   title      = {Solving Multi-task Optimization Problems with Adaptive Knowledge Transfer via Anomaly Detection},
%   author     = {Wang, Chao and Liu, Jing and Wu, Kai and Wu, Zhaoyang},
%   journal    = {IEEE Transactions on Evolutionary Computation},
%   year       = {2021},
%   pages      = {1-1},
%   doi        = {10.1109/TEVC.2021.3068157},
% }
%--------------------------------------------------------------------------

%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 Yanchi Li. You are free to use the MTO-Platform for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "MTO-Platform" and cite
% or footnote "https://github.com/intLyc/MTO-Platform"
%--------------------------------------------------------------------------

properties (SetAccess = private)
    TRP = 0.1
    MuC = 2
    MuM = 5
    Alpha=1
end

methods
    function Parameter = getParameter(Algo)
        Parameter = {'TRP: Probability of the Knowledge Transfer', num2str(Algo.TRP), ...
                'MuC: Simulated Binary Crossover', num2str(Algo.MuC), ...
                'MuM: Polynomial Mutation', num2str(Algo.MuM),...
                'Alpha',num2str(Algo.Alpha)};
    end

    function Algo = setParameter(Algo, Parameter)
        i = 1;
        Algo.TRP = str2double(Parameter{i}); i = i + 1;
        Algo.MuC = str2double(Parameter{i}); i = i + 1;
        Algo.MuM = str2double(Parameter{i}); i = i + 1;
        Algo.Alpha = str2double(Parameter{i}); i = i + 1;
    end

    function run(Algo, Prob)
        % Initialization
        population = Initialization(Algo, Prob, Individual);
        epsilon = zeros(1, Prob.T); % Parameter of the anomaly detection model

        alpha = [];
        allmodels = {};
        for t = 1:Prob.T
            % initialize allmodals
            model = ProbabilityModel('mvarnorm'); 
            model = ProbabilityModel.buildmodel(model,population{t}.Decs);
            allmodels{t} = model;
        end

        while Algo.notTerminated(Prob)
             if Algo.FE <= Prob.maxFE * Algo.Alpha
            for t = 1:Prob.T
                offs = Individual.empty();
                mmodel = MixModel({allmodels{1:t-1},allmodels{t+1:end}}); 
                mmodel = MixModel.createtable(mmodel, population{t}.Decs, true, 'mvarnorm'); 
                mmodel = MixModel.EMstacking(mmodel); 
                mmodel = MixModel.mutate(mmodel); 
                alpha = [alpha; mmodel.alpha];
                soffs = MixModel.sample(mmodel, length(population{t})); 
                for i = 1:size(soffs,1)
                    offs(i).Dec = soffs(i,:);
                    offs(i).Dec(offs(i).Dec > 1) = 1;
                    offs(i).Dec(offs(i).Dec < 0) = 0;
                end
                offs = Algo.Evaluation(offs, Prob, t);
                population{t} = Selection_Tournament(population{t}, offs);
            end
             end
            for t = 1:Prob.T
                % Generation
                offspring = Algo.Generation(population{t});

                % Knowledge Transfer
                if rand() < Algo.TRP
                    if Algo.Gen == 1
                        NL = 1;
                    else
                        NL = epsilon(t);
                    end

                    kpool = 1:Prob.T; kpool(kpool == t) = [];
                    kpool = kpool(randperm(Prob.T - 1, min(Prob.T - 1, 10))); % Reduce time complexity in many-task
                    his_pop_dec = [];
                    for k = kpool
                        his_pop_dec = [his_pop_dec; population{k}.Decs];
                    end

                    tfsol = Algo.learn_anomaly_detection(offspring.Decs, his_pop_dec, NL);

                    transfer_pop = Individual.empty();
                    for i = 1:size(tfsol, 1)
                        c = Individual();
                        c.Dec = tfsol(i, :);
                        c.Dec(c.Dec > 1) = 1;
                        c.Dec(c.Dec < 0) = 0;
                        transfer_pop = [transfer_pop, c];
                    end

                    % Evaluation
                    offspring = Algo.Evaluation(offspring, Prob, t);
                    transfer_pop = Algo.Evaluation(transfer_pop, Prob, t);
                    % Selection
                    [population{t}, rank] = Selection_Elit(population{t}, [offspring, transfer_pop]);
                    succ_num = sum(rank(1:length(population{t})) > length(population{t}) + length(offspring));
                    % Parameter adaptation strategy via elitism
                    epsilon(t) = succ_num ./ size(tfsol, 1);
                else
                    % Evaluation
                    offspring = Algo.Evaluation(offspring, Prob, t);
                    % Selection
                    population{t} = Selection_Elit(population{t}, offspring);
                end
            end
        end
    end

    function offspring = Generation(Algo, population)
        indorder = randperm(length(population));
        count = 1;
        for i = 1:ceil(length(population) / 2)
            p1 = indorder(i);
            p2 = indorder(i + fix(length(population) / 2));
            offspring(count) = population(p1);
            offspring(count + 1) = population(p2);

            [offspring(count).Dec, offspring(count + 1).Dec] = GA_Crossover(population(p1).Dec, population(p2).Dec, Algo.MuC);

            offspring(count).Dec = GA_Mutation(offspring(count).Dec, Algo.MuM);
            offspring(count + 1).Dec = GA_Mutation(offspring(count + 1).Dec, Algo.MuM);

            for x = count:count + 1
                offspring(x).Dec(offspring(x).Dec > 1) = 1;
                offspring(x).Dec(offspring(x).Dec < 0) = 0;
            end
            count = count + 2;
        end
    end

    function tfsol = learn_anomaly_detection(Algo, curr_pop, his_pop, NL)
        %% Learning anomaly detection model of task tn
        % Input: curr_pop (Dec matrix), his_pop (Dec matrix), NL (anomaly detection parameter)
        % Output: tfsol (candidate transferred solutions)

        % Sample, make sure that the fitted covariance is a square, symmetric, positive definite matrix.
        nsamples = floor(0.01 * size(curr_pop, 1));
        randMat = rand(nsamples, size(curr_pop, 2));
        curr_pop = [curr_pop; randMat];

        % Fit
        mmean = mean(curr_pop);
        sstd1 = cov(curr_pop);
        sstd = sstd1 + (10e-6) * eye(size(curr_pop, 2));

        % Calculate the scores
        [Dec, ~] = unique(his_pop, 'rows');
        Y = mvnpdf(Dec(:, 1:size(curr_pop, 2)), mmean, sstd);

        % Select the candidate transferred solutions
        [~, ii] = sort(Y, 'descend');
        if NL == 0
            mm = Y(1); % Ensure that the number of transferred individuals is not 0
        else
            mm = Y(ii(ceil(size(Y, 1) * NL)));
        end

        % Count the number of candidate transferred solutions
        tte = Y >= mm;
        tfsol = Dec(tte, :);
    end
end
end
