classdef ProbabilityModel % Works reliably for 2(+) Dimensional distributions
    properties
        modeltype; % multivariate normal ('mvarnorm' - for real coded) or univariate marginal distribution ('umd' - for binary coded)
        %连续
        mean_noisy;
        mean_true;
        covarmat_noisy;
        covarmat_true;
        %离散 二元编码
        probofone_noisy;
        probofone_true; %为1的概率
        probofzero_noisy;
        probofzero_true;    %为0的概率
        vars;
    end
    methods (Static)
        function model = ProbabilityModel(type)
            model.modeltype = type;
        end
        function solutions = sample(model,nos)
            if strcmp(model.modeltype,'mvarnorm')
                %solutions = mvnrnd(model.mean_true,model.covarmat_true,nos);
                solutions = mvnrnd(model.mean_noisy,model.covarmat_noisy,nos);
            elseif strcmp(model.modeltype,'umd')
                solutions = rand(nos,model.vars);
                for i = 1:nos
                    index1 = solutions(i,:) <= model.probofone_true;
                    index0 = solutions(i,:) > model.probofone_true;
                    solutions(i,index1) = 1;
                    solutions(i,index0) = 0;
                end
            end
        end
        function probofsols = pdfeval(model,solutions)
            %% 根据模型计算某个解个体的概率
            if strcmp(model.modeltype,'mvarnorm')
                probofsols = mvnpdf(solutions,model.mean_noisy,model.covarmat_noisy);
            elseif strcmp(model.modeltype,'umd')
                nos = size(solutions,1);
                probofsols = zeros(nos,1);
                probvector = zeros(1,model.vars);
                for i = 1:nos
                    index = solutions(i,:) == 1;
                    probvector(index) = model.probofone_noisy(index);
                    index = solutions(i,:) == 0;
                    probvector(index) = model.probofzero_noisy(index);
                    probofsols(i) = prod(probvector);%返回各列元素概率的乘积
                end
            end
        end

        %% 根据传入的解种群构建模型
        function model = buildmodel(model,solutions)
            [pop,model.vars] = size(solutions);
            if strcmp(model.modeltype,'mvarnorm')
                model.mean_true = mean(solutions);
                covariance = cov(solutions);
                model.covarmat_true = diag(diag(covariance)); % Simplifying to univariate distribution by ignoring off diagonal terms of covariance matrix

                solutions_noisy = [solutions;rand(round(0.1*pop),model.vars)];
                model.mean_noisy = mean(solutions_noisy);
                covariance = cov(solutions_noisy);
                model.covarmat_noisy = diag(diag(covariance));% Simplifying to univariate distribution by ignoring off diagonal terms of covariance matrix
                model.covarmat_noisy = cov(solutions_noisy);
            elseif strcmp(model.modeltype,'umd')
                model.probofone_true = mean(solutions);%表示每列元素为1的真实概率，返回数组每列均值的行向量
                model.probofzero_true = 1 - model.probofone_true;%每列元素为0的真实概率
                solutions_noisy = [solutions;round(rand(round(0.1*pop),model.vars))];%添加噪声
                model.probofone_noisy = mean(solutions_noisy);%添加噪声后每列元素为1的概率
                model.probofzero_noisy = 1 - model.probofone_noisy;%添加噪声后每列元素为0的概率
            end
        end
    end
end