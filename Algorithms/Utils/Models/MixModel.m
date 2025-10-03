classdef MixModel % Works reliably for 2(+) Dimensional distributions
    properties
        model_list; % cell array of ProbabilityModels
        alpha; % weights of the models in stacking of mixture models
        noms; % number of models
        probtable; % Probability table required for stacking EM algorithm
        nsols; % number of solutions in probability table
    end
    methods (Static)
        function mmodel = MixModel(allmodels)
            %% 初始化混合模型
            mmodel.model_list = allmodels; %模型列表
            mmodel.noms = length(allmodels);%子模型数量
            mmodel.alpha = (1/mmodel.noms)*ones(1,mmodel.noms);%各子模型的概率
        end
        function mmodel = EMstacking(mmodel,varargin)
            
            % 根据所得概率表，进行EM迭代过程，使参数趋近真实概率分布的参数
            n = numel(varargin);
            if n == 0
                 iterations = 30;
            elseif n == 1
                iterations = varargin{1};
            else
                return;
            end

            for i = 1:iterations
                talpha = mmodel.alpha;
                probvector = mmodel.probtable*talpha';
                for j = 1:mmodel.noms
                    talpha(j) = sum((1/mmodel.nsols)*talpha(j)*mmodel.probtable(:,j)./probvector);
                end
                mmodel.alpha = talpha;
            end
        end
        function mmodel = mutate(mmodel)
            modifalpha = max(mmodel.alpha+normrnd(0,0.01,[1,mmodel.noms]),0); %%%%%%%% Determining std dev for mutation can be a parameteric study %%%%%%%%%%%%%%%%
            pusum = sum(modifalpha);
            if pusum == 0 % Then complete weightage assigned to target model alone
                mmodel.alpha = zeros(1,mmodel.noms);
                mmodel.alpha(mmodel.noms) = 1;
            else
                mmodel.alpha = modifalpha/pusum;
            end
        end
        function solutions = sample(mmodel,nos)
            indsamples = ceil(nos*mmodel.alpha);%得到总共取样多少个体
            totalsamples = sum(indsamples);
            solutions = [];
            for i = 1:mmodel.noms
                if indsamples(i) == 0
                    continue;
                else
                    sols = ProbabilityModel.sample(mmodel.model_list{i},indsamples(i));
                    solutions = [solutions; sols];
                end
            end
            solutions = solutions(randperm(totalsamples),:);
            solutions = solutions(1:nos,:);
        end

        function mmodel = createtable(mmodel,solutions,CV,type)
            %% 建表
            if CV
                mmodel.noms = mmodel.noms+1; %%%%%% NOTE: Last model in the list is the target model
                mmodel.model_list{mmodel.noms} = ProbabilityModel(type); %初始化一个目标模型
                mmodel.model_list{mmodel.noms} = ProbabilityModel.buildmodel(mmodel.model_list{mmodel.noms},solutions);%利用目标任务种群训练目标模型
                mmodel.alpha = (1/mmodel.noms)*ones(1,mmodel.noms);%各子模型的初始概率
                nos = size(solutions,1);%解的数量nos
                mmodel.probtable = ones(nos,mmodel.noms);%初始化一个表  解的个数*模型数量
                for j =1:mmodel.noms-1
                    % 表中前两列元素为 将种群个体代入模型后得到的概率？？？？
                    % 将已有样本（目标任务种群）代入初始的模型，求出发生该情况（个体为目标任务中群）的概率
                    % 遍历，得到概率表
                    mmodel.probtable(:,j) = ProbabilityModel.pdfeval(mmodel.model_list{j},solutions);
                end
                for i = 1:nos % 保留一个个体作为交叉验证方案   Leave-one-out cross validation scheme 
                    x = [solutions(1:i-1,:);solutions(i+1:nos,:)]; %利用其他行构建模型
                    tmodel = ProbabilityModel(type);
                    tmodel = ProbabilityModel.buildmodel(tmodel,x);
                    %将保留个体代入模型得到表中第三列概率
                    mmodel.probtable(i,mmodel.noms) = ProbabilityModel.pdfeval(tmodel,solutions(i,:)); %剩余那行用于计算概率
                end
            else
                nos = size(solutions,1);
                mmodel.probtable = ones(nos,mmodel.noms);
                for j =1:mmodel.noms
                    mmodel.probtable(:,j) = ProbabilityModel.pdfeval(mmodel.model_list{j},solutions);
                end
            end
            mmodel.nsols = nos;
        end
    end
end