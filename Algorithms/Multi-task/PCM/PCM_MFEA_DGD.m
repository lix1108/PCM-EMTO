classdef PCM_MFEA_DGD < Algorithm
% <Multi-task> <Single-objective> <None>

%------------------------------- Reference --------------------------------
% @Article{Liu2023MFEA-DGD,
%   author   = {Liu, Zhaobo and Li, Guo and Zhang, Haili and Liang, Zhengping and Zhu, Zexuan},
%   journal  = {IEEE Transactions on Cybernetics},
%   title    = {Multifactorial Evolutionary Algorithm Based on Diffusion Gradient Descent},
%   year     = {2023},
%   pages    = {1-13},
% }
%--------------------------------------------------------------------------

%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 Yanchi Li. You are free to use the MTO-Platform for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "MTO-Platform" and cite
% or footnote "https://github.com/intLyc/MTO-Platform"
%--------------------------------------------------------------------------

properties (SetAccess = private)
    RMP = 0.3
    Gamma = 0.1
    Alpha = 0.25
end

methods
    function Parameter = getParameter(Algo)
        Parameter = {'RMP: Random Mating Probability', num2str(Algo.RMP), ...
                'Gamma', num2str(Algo.Gamma),...
                'Alpha', num2str(Algo.Alpha)};
    end

    function setParameter(Algo, Parameter)
        i = 1;
        Algo.RMP = str2double(Parameter{i}); i = i + 1;
        Algo.Gamma = str2double(Parameter{i}); i = i + 1;
        Algo.Alpha = str2double(Parameter{i}); i = i + 1;
    end

    function run(Algo, Prob)
        % Initialize
        population = Initialization_MF(Algo, Prob, Individual_MF);
        sigma = 0;

        alpha = [];
        allmodels = {};
        all_offs = Individual_MF.empty();

        for t = 1:Prob.T
            % initialize Models
            model = ProbabilityModel('mvarnorm'); 
            model = ProbabilityModel.buildmodel(model,population([population.MFFactor] == t).Decs);
            allmodels{t} = model;
        end

        while Algo.notTerminated(Prob)
            for t = 1:Prob.T
                if Algo.FE <= Prob.maxFE * Algo.Alpha
                offs = Individual_MF.empty();
                mmodel = MixModel({allmodels{1:t-1},allmodels{t+1:end}}); 
                mmodel = MixModel.createtable(mmodel, population([population.MFFactor] == t).Decs, true, 'mvarnorm'); 
                mmodel = MixModel.EMstacking(mmodel); 
                mmodel = MixModel.mutate(mmodel); 
                alpha = [alpha; mmodel.alpha];
                soffs = MixModel.sample(mmodel, length(population([population.MFFactor] == t))); 
                for i = 1:size(soffs,1)
                    offs(i).Dec = soffs(i,:);
                    offs(i).Dec(offs(i).Dec > 1) = 1;
                    offs(i).Dec(offs(i).Dec < 0) = 0;
                    offs(i).MFFactor = t;
                end
                offs = Algo.Evaluation(offs, Prob, t);
                for i = 1:length(offs)
                    offs(i).MFObj = inf(1, Prob.T);
                    offs(i).MFCV = inf(1, Prob.T);
                    offs(i).MFObj(t) = offs(i).Obj;
                    offs(i).MFCV(t) = offs(i).CV;
                end
                all_offs = [all_offs, offs];
            end 
            % Selection
            population = Selection_MF(population, all_offs, Prob);        
            end
            f = randperm(5);
            if sigma <= 0
                for i = 1:5
                    if f(1) == i
                        sigma = 10^(-i);
                    end
                end
            end
            for t = 1:Prob.T
                Dect = population([population.MFFactor] == t).Decs;
                maxDec{t} = max(Dect);
                minDec{t} = min(Dect);
            end
            % Generation
            [offspring, offspring2] = Algo.Generation(Prob, population, sigma, maxDec, minDec);
            % Evaluation
            offspring_temp = Individual_MF.empty();
            for t = 1:Prob.T
                offspring_t = offspring([offspring.MFFactor] == t);
                offspring_t = Algo.Evaluation(offspring_t, Prob, t);
                offspring_t = [offspring_t, offspring2([offspring2.MFFactor] == t)];
                for i = 1:length(offspring_t)
                    offspring_t(i).MFObj = inf(1, Prob.T);
                    offspring_t(i).MFCV = inf(1, Prob.T);
                    offspring_t(i).MFObj(t) = offspring_t(i).Obj;
                    offspring_t(i).MFCV(t) = offspring_t(i).CV;
                end
                offspring_temp = [offspring_temp, offspring_t];
            end
            offspring = offspring_temp;
            % Selection
            population = Selection_MF(population, offspring, Prob);
        end
    end

    function [offspring, offspring2] = Generation(Algo, Prob, population, sigma, maxDec, minDec)
        L = 0;
        indorder = randperm(length(population));
        count = 1;
        for i = 1:ceil(length(population) / 2)
            p1 = indorder(i);
            p2 = indorder(i + fix(length(population) / 2));
            offspring(count) = population(p1);
            offspring(count + 1) = population(p2);
            offspring2(count) = population(p1);
            offspring2(count + 1) = population(p2);
            offspring2(count + 2) = population(p1);
            offspring2(count + 3) = population(p2);

            k = 0.7 + 0.6 * rand();
            QWE = zeros(2, length(population(p1).Dec));
            E = RandOrthMat(length(population(p1).Dec), 1);
            p = [p1, p2];
            factor = [population(p).MFFactor];
            for x = 1:2
                sd = E';
                offspring2(count + 2 * (x - 1)).Dec = population(p(x)).Dec + sd(1, :) .* sigma;
                offspring2(count + 2 * (x - 1) + 1).Dec = population(p(x)).Dec - sd(1, :) .* sigma;
                offspring2(count + 2 * (x - 1)) = Algo.Evaluation(offspring2(count + 2 * (x - 1)), Prob, factor(x));
                offspring2(count + 2 * (x - 1) + 1) = Algo.Evaluation(offspring2(count + 2 * (x - 1) + 1), Prob, factor(x));
                L1 = offspring2(count + 2 * (x - 1)).Obj - offspring2(count + 2 * (x - 1) + 1).Obj;
                QWE(x, :) = QWE(x, :) + (sd .* L1) / sigma;
            end

            if norm(QWE) > L
                L = (1 - Algo.Gamma) * norm(QWE) + Algo.Gamma * L;
            end

            if (population(p1).MFFactor == population(p2).MFFactor) || rand() < Algo.RMP
                % crossover
                p = [p1, p2]; r1 = randi(2); r2 = mod(r1, 2) + 1;
                t = [population(p).MFFactor];
                factor = t(randi(2));
                offspring(count).Dec = Algo.Crossover(population(p(r1)).Dec, population(p(r2)).Dec, QWE, L, sigma);
                if rand() > mod(Algo.Gen, 2) % OBL
                    offspring(count + 1).Dec = 1 - offspring(count).Dec;
                else
                    offspring(count + 1).Dec = k * (maxDec{factor} + minDec{factor}) - offspring(count).Dec;

                end
                % imitation
                offspring(count).MFFactor = population(p(randi(2))).MFFactor;
                offspring(count + 1).MFFactor = population(p(randi(2))).MFFactor;
            else
                % mutation
                offspring(count).Dec = Algo.Mutation(population(p1).Dec, QWE(1, :), L, sigma);
                offspring(count + 1).Dec = Algo.Mutation(population(p2).Dec, QWE(2, :), L, sigma);
                % imitation
                offspring(count).MFFactor = population(p1).MFFactor;
                offspring(count + 1).MFFactor = population(p2).MFFactor;
            end
            for x = count:count + 1
                offspring(x).Dec(offspring(x).Dec > 1) = 1;
                offspring(x).Dec(offspring(x).Dec < 0) = 0;
            end
            count = count + 2;
        end
    end

    function [OffDec1, OffDec2] = Crossover(Algo, ParDec1, ParDec2, QWE, L, sigma)
        D = length(ParDec1);
        u = rand(1, D);
        cf = zeros(1, D);
        cf(u <= 0.5) = 0.6 * rand();
        cf(u > 0.5) = -0.6 * rand();

        ParDec1 = ParDec1 - QWE(1, :) .* sigma / L;
        ParDec2 = ParDec2 - QWE(2, :) .* sigma / L;
        OffDec1 = 0.5 * ((1 + cf) .* ParDec1 + (1 - cf) .* ParDec2);
        OffDec2 = 0.5 * ((1 + cf) .* ParDec2 + (1 - cf) .* ParDec1);
    end

    function Dec = Mutation(Algo, Dec, qq, L, sigma)
        Dec = Dec - qq .* sigma / L;
    end
end
end
