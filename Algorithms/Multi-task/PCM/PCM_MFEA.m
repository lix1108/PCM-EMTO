classdef PCM_MFEA < Algorithm
% <Multi-task> <Single-objective> <None/Constrained>

%------------------------------- Reference --------------------------------
% @Article{Gupta2016MFEA,
%   title      = {Multifactorial Evolution: Toward Evolutionary Multitasking},
%   author     = {Gupta, Abhishek and Ong, Yew-Soon and Feng, Liang},
%   journal    = {IEEE Transactions on Evolutionary Computation},
%   year       = {2016},
%   number     = {3},
%   pages      = {343-357},
%   volume     = {20},
%   doi        = {10.1109/TEVC.2015.2458037},
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
    MuC = 2
    MuM = 5
    iter = 30
    Alpha = 0.25
end

methods
    function Parameter = getParameter(Algo)
        Parameter = {'RMP: Random Mating Probability', num2str(Algo.RMP), ...
                'MuC: Simulated Binary Crossover', num2str(Algo.MuC), ...
                'MuM: Polynomial Mutation', num2str(Algo.MuM), ...
                'iter: EM iteration', num2str(Algo.iter),...
                'Alpha', num2str(Algo.Alpha)};
    end

    function setParameter(Algo, Parameter)
        i = 1;
        Algo.RMP = str2double(Parameter{i}); i = i + 1;
        Algo.MuC = str2double(Parameter{i}); i = i + 1;
        Algo.MuM = str2double(Parameter{i}); i = i + 1;
        Algo.iter = str2double(Parameter{i}); i = i + 1;
        Algo.Alpha = str2double(Parameter{i}); i = i + 1;
    end

    function run(Algo, Prob)
        % Initialize
        population = Initialization_MF(Algo, Prob, Individual_MF);

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
                    mmodel = MixModel.EMstacking(mmodel,Algo.iter); 
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
                % Selection
                population = Selection_MF(population, all_offs, Prob);
                end
            end
            % Generation
            offspring = Algo.Generation(population);
            % Evaluation
            offspring_temp = Individual_MF.empty();
            for t = 1:Prob.T
                offspring_t = offspring([offspring.MFFactor] == t);
                offspring_t = Algo.Evaluation(offspring_t, Prob, t);
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

    function offspring = Generation(Algo, population)
        indorder = randperm(length(population));
        count = 1;
        for i = 1:ceil(length(population) / 2)
            p1 = indorder(i);
            p2 = indorder(i + fix(length(population) / 2));
            offspring(count) = population(p1);
            offspring(count + 1) = population(p2);

            if (population(p1).MFFactor == population(p2).MFFactor) || rand() < Algo.RMP
                % crossover
                [offspring(count).Dec, offspring(count + 1).Dec] = GA_Crossover(population(p1).Dec, population(p2).Dec, Algo.MuC);
                % imitation
                p = [p1, p2];
                offspring(count).MFFactor = population(p(randi(2))).MFFactor;
                offspring(count + 1).MFFactor = population(p(randi(2))).MFFactor;
            else
                % mutation
                offspring(count).Dec = GA_Mutation(population(p1).Dec, Algo.MuM);
                offspring(count + 1).Dec = GA_Mutation(population(p2).Dec, Algo.MuM);
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
end
end
