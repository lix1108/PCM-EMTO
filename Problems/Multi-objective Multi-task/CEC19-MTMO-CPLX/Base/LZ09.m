classdef LZ09
%UNTITLED4 �˴���ʾ�йش����ժҪ
%   �˴���ʾ��ϸ˵��

properties
    dim; %���߱���ά��
    numOfObjective; %Ŀ����
    ltype;
    dtype;
    ptype;
end

methods
    %% functionname: function description
    function obj = LZ09(dim, numOfObjective, ltype, dtype, ptype)
        obj.dim = dim;
        obj.numOfObjective = numOfObjective;
        obj.ltype = ltype;
        obj.dtype = dtype;
        obj.ptype = ptype;
    end

    function [beta] = psfunc(obj, x, x1, css) %����type ��ltype
        %x  : �ֺ���ż�ľ��߱���    x1:���߱����ĵ�һά
        %obj.dim������ά��   obj.ltype���������� ��css��class of index ;  xective_num��Ŀ����
        if (obj.ltype == 21)
            if mod(size(x, 2), 2) == 0
                qq = [3:2:obj.dim];
            else
                qq = [2:2:obj.dim];
            end
            x = 2 .* (x - 0.5);
            beta = x - x1.^(0.5 .* (obj.dim + 3 .* qq - 8) ./ (obj.dim - 2));
        end

        if (obj.ltype == 22)
            if mod(size(x, 2), 2) == 0
                qq = [3:2:obj.dim];
            else
                qq = [2:2:obj.dim];
            end
            x = 2 .* (x - 0.5);
            theta = sin(6 * pi * x1 + (pi .* qq) ./ obj.dim);
            beta = x - theta;
        end

        if obj.ltype == 23

            if mod(size(x, 2), 2) == 0
                qq = [3:2:obj.dim];
            else
                qq = [2:2:obj.dim];
            end
            theta = 6 * pi * x1 + (pi .* qq) ./ obj.dim;
            ra = 0.8 * x1;
            x = 2 .* (x - 0.5);
            if (css == 1)

                beta = x - ra .* cos(theta);
            else
                beta = x - ra .* sin(theta);
            end
        end

        if obj.ltype == 24
            if mod(size(x, 2), 2) == 0
                qq = [3:2:obj.dim];
            else
                qq = [2:2:obj.dim];
            end
            theta = 6 * pi * x1 + (pi .* qq) ./ obj.dim;
            ra = 0.8 * x1;
            x = 2 .* (x - 0.5);
            if (css == 1)
                beta = x - ra .* cos(theta ./ 3);
            else
                beta = x - ra .* sin(theta);
            end
        end

        if obj.ltype == 25
            rho = 0.8;
            phi = pi * x1;
            if mod(size(x, 2), 2) == 0
                qq = [3:2:obj.dim];
            else
                qq = [2:2:obj.dim];
            end
            theta = 6 * pi * x1 + (pi .* qq) ./ obj.dim;
            x = 2 .* (x - 0.5);
            if css == 1
                beta = x - rho * sin(phi) .* sin(theta);
            elseif css == 2
                beta = x - rho * sin(phi) .* cos(theta);
            else
                beta = x - rho * cos(phi);
            end
        end

        if obj.ltype == 26
            if mod(size(x, 2), 2) == 0
                qq = [3:2:obj.dim];
            else
                qq = [2:2:obj.dim];
            end
            theta = 6 * pi * x1 + (pi .* qq) ./ obj.dim;
            ra = 0.3 * x1 .* (x1 .* cos(4 .* theta) + 2);
            x = 2 .* (x - 0.5);
            if (css == 1)
                beta = x - ra .* cos(theta);
            else
                beta = x - ra .* sin(theta);
            end
        end
    end
    %-----------------psfun  end-----------------------

    %-----------------psfun3 begin -----------------------
    %% psfunc3: function description
    function beta = psfunc3(obj, x, x1, x2, order)
        % obj :����ʵ��		x��J1 , J2 , J3 �� order ��J1 , J2, J3��ԭ���߱������е����
        if obj.ltype == 32
            theta = 2 * pi * x1 + pi .* order ./ obj.dim;
            x = 4 .* (x - 0.5);
            beta = x - 2 * x2 .* sin(theta);
        end
    end
    % ---------- psfun3 end---------------------------------

    function alpha = alphaFunction(obj, x) %�Ӻ�����������һ����

        %inpu:obj,�����ʵ��������  ,x:����ľ��߱���   ����type����ptype
        %output ��alpha ����Ŀ�꺯����ͬ������㲿��
        if obj.numOfObjective == 2 %��Ŀ������
            switch obj.ptype
                case 21
                    alpha(:, 1) = x(:, 1);
                    alpha(:, 2) = 1 - sqrt(x(:, 1));
                case 22
                    alpha(:, 1) = x(:, 1);
                    alpha(:, 2) = 1 - x(:, 1).^2;
                case 23
                    alpha(:, 1) = x(:, 1);
                    alpha(:, 2) = 1 - sqrt(x(:, 1)) - x(:, 1) * sin(10 * x(:, 1).^2 * pi);
                case 24
                    alpha(:, 1) = x(:, 1);
                    alpha(:, 2) = 1 - x(:, 1) - 0.05 * sin(4 * pi * x(:, 1));
            end
        else %��Ŀ������
            switch obj.ptype
                case 31
                    alpha(:, 1) = cos(x(:, 1) * (pi / 2)) .* cos(x(:, 2) * (pi / 2));
                    alpha(:, 2) = cos(x(:, 1) * (pi / 2)) .* sin(x(:, 2) * (pi / 2));
                    alpha(:, 3) = sin(x(:, 1) * (pi / 2));
                case 32
                    alpha(:, 1) = 1 - cos(x(:, 1) * (pi / 2)) .* cos(x(:, 2) * (pi / 2));
                    alpha(:, 2) = 1 - cos(x(:, 1) * (pi / 2)) .* sin(x(:, 2) * (pi / 2));
                    alpha(:, 3) = 1 - sin(x(:, 1) * (pi / 2));
                case 33
                    alpha(:, 1) = x(:, 1);
                    alpha(:, 2) = x(:, 2);
                    alpha(:, 3) = 3 - (sin(3 * pi * x(:, 1)) + sin(3 * pi * x(:, 2))) - 2 * (x(:, 1) + x(:, 2));
                case 34
                    alpha(:, 1) = x(:, 1) .* x(:, 2);
                    alpha(:, 2) = x(:, 1) .* (1 - x(:, 1));
                    alpha(:, 3) = 1 - x(:, 1);
            end
        end
    end
    % --------------------------------------------------------
    function beta = betaFunction(oddEven_x, obj)
        %input: oddEven_x : ��ʱֻ�Ǵ�������ά����ż��ά�ľ��߱�������������ά��������һά
        %output : beta����������ۼӺͲ��֣�ֻ���ۼӵĲ��ֲ�ͬ obj.dtype��ֵ��ͬ
        beta = zeros(size(oddEven_x, 1));
        dim = size(oddEven_x, 2);
        switch obj.dtype
            case 1
                beta = sum(oddEven_x.^2, 2);
                beta = 2 * beta / dim;
            case 2
                a = [1:1:dim];
                a = sqrt (a);
                oddEven_x = oddEven_x.^2;
                beta = sum(oddEven_x .* a, 2);
                beta = 2 * beta / dim;
            case 3
                oddEven_x1 = oddEven_x; %������ҪoddEven_x��һ������
                beta = sum((2 .* oddEven_x).^2, 2) - sum(cos(4 * pi * oddEven_x1), 2);
                beta = beta + dim;
                beta = 2 * beta / dim;
            case 4
                sum1 = 0;
                prod = 1;
                oddEven_x1 = oddEven_x; %oddEven_x :�����ۼ�   oddEven_x1�������۳�
                sum1 = sum((2 .* oddEven_x).^2, 2);
                a = [1:dim];
                a = sqrt(a);
                prod = cos(10 * pi * (2 .* oddEven_x1) ./ a);
                prod = cumprod(prod);
                a = prod(end);
                beta = 2 * (sum1 - 2 * a + 2) / dim;
        end

    end
    %---------------------------betaFunction finished--------------------------

    %---------------------------objectiveFunction begin-----------------------
    %% objectiveFunction:nction description
    function fitness = objectiveFunction(obj, x)
        % input: x , ���߱���  obj������ʵ��
        % output: fitness_1 :Ŀ�꺯�� 1������ֵ �� Ŀ�꺯��2������ֵ
        ltypeTable = [21, 22, 23, 24, 26];
        if obj.numOfObjective == 2 %��Ŀ������
            if ismember(obj.ltype, ltypeTable)
                J1 = x(:, 3:2:end);
                a = psfunc(obj, J1, x(:, 1), 1);
                %��ȡ��ż��ά
                J2 = x(:, 2:2:end);
                b = psfunc(obj, J2, x(:, 1), 2);
                g = betaFunction(a, obj);
                h = betaFunction(b, obj);
                alpha1 = alphaFunction(obj, x);
                fitness_1 = alpha1(:, 1) + h;
                fitness_2 = alpha1(:, 2) + g;
                fitness = [fitness_1, fitness_2];
            else %ltype ��25���������� %�Ѿ��߱�����Ϊ3��
                J1 = x(:, 4:3:end);
                a = psfunc(obj, J1, x(:, 1), 1);
                J2 = x(:, 2:3:end);
                b = psfunc(obj, J2, x(:, 1), 2);
                J3 = x(:, 3:3:end);
                c = psfunc(obj, J3, x(:, 1), 3);
                % �������Ƹ�A ,ż�����Ƹ�B
                a = [a, c(:, 1:2:end)];
                b = [b, b(:, 2:2:end)];
                g = betaFunction(a, obj);
                h = betaFunction(b, obj);
                alpha1 = alphaFunction(obj, x);
                fitness_1 = alpha1(:, 1) + h;
                fitness_2 = alpha1(:, 2) + g;
                fitness = [fitness_1, fitness_2];
            end
        else

            J1 = x(:, 4:3:end);
            J2 = x(:, 5:3:end);
            J3 = x(:, 3:3:end);
            order1 = [4:3:size(x, 2)];
            order2 = [5:3:size(x, 2)];
            order3 = [3:3:size(x, 2)];
            a = psfunc3(obj, J1, x(:, 1), x(:, 2), order1);
            b = psfunc3(obj, J2, x(:, 1), x(:, 2), order2);
            c = psfunc3(obj, J3, x(:, 1), x(:, 2), order3);
            g = betaFunction(a, obj);
            h = betaFunction(b, obj);
            e = betaFunction(c, obj);
            alpha1 = alphaFunction(obj, x);
            fitness_1 = alpha1(:, 1) + g;
            fitness_2 = alpha1(:, 2) + h;
            fitness_3 = alpha1(:, 3) + e;
            fitness = [fitness_1, fitness_2, fitness_3];
        end
    end
end
end
