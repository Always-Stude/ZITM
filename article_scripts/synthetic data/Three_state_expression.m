clc
clear all
close all

% 基础参数
N1=500;
GN=10000;
delta=1;
T = 200;
t5=T;

numWorkers = 10;  % 对应JSUB -n 10
% 初始化存储矩阵
parameter = cell(GN+1, 8);
expression = cell(GN+1, N1+1);  
parameter(1, :) = {'', 'lam1', 'lam2', 'gamma', 'nu', '', 'bf', 'bs'};
expression{1,1} = '';
for i=1:N1
    expression{1, i+1} = sprintf('cell%d', i);
end
for i=1:GN
    parameter{i+1,1}=sprintf('gene%d', i);
    expression{i+1,1}=sprintf('gene%d', i);
end
parameter(2:end, 2:end) = num2cell(NaN(GN, 7));
expression(2:end, 2:end) = num2cell(NaN(GN, N1));

% 预分配并行临时存储
par_lam1 = zeros(GN,1);
par_lam2 = zeros(GN,1);
par_gamma = zeros(GN,1);
par_nu = zeros(GN,1);
par_bf = zeros(GN,1);
par_bs = zeros(GN,1);
par_S5 = zeros(GN,N1);

% 开启本地并行池，内部并行循环基因
parpool('local');
parfor zushu = 1:GN
    rng(zushu); %   固定随机种子，结果可复现
    % 随机参数
    lam = 0.2 + 5.8 * rand();
    lam1 = lam;
    lam2 = lam;  
    gamma = 0.1 + 2.9 * rand(); 
    nu = 10 + 20 * rand();
    nu1=nu;
    
    S5 = [];
    % 单基因多个样本模拟
    for i = 1:N1
        X(1) = 0; t(1) = 0;
        s(1)=1;
        n = 1; 
        Xall = [X(1)]; Sall = [s(1)]; tall = [t(1)];
        
        while t(n) <= T
            a1 = nu1;
            a2 = gamma;
            a3 = lam1;
            a4 = lam2;
            a5 = X(n)*delta;
            n = n+1;
            
            if s(n-1) == 0
                a0 = a1+a2+a5;
                r1 = rand; r2=rand;
                tau = -log(r1)/a0;
                if a0*r2 <= a1
                    X(n)=X(n-1)+1;  s(n)=0;
                elseif a0*r2 <=a1+a2
                    X(n)=X(n-1);  s(n)=1;
                else
                    X(n)=X(n-1)-1; s(n)=0;
                end
            elseif s(n-1) == 1
                a0 = a3+a5;
                r1=rand;  r2=rand;
                tau=-log(r1)/a0;
                if a0*r2 <= a3
                    X(n)=X(n-1); s(n)=2;
                else
                    X(n)=X(n-1)-1; s(n)=1;
                end
            else
                a0 = a4+a5;
                r1=rand;  r2=rand;
                tau=-log(r1)/a0;
                if a0*r2 <= a4
                    X(n)=X(n-1); s(n)=0;
                else
                    X(n)=X(n-1)-1; s(n)=2;
                end
            end
            t(n) = t(n-1) + tau;
            Xall = [Xall X(n)]; Sall = [Sall s(n)]; tall = [tall t(n)];
        end
        if t(n)<T
            Xall=[Xall X(n)]; Sall=[Sall s(n)]; tall=[tall T];
        end
        % 截取稳态数值
        val = NaN;
        for j=2:n
            if tall(j-1) <=t5 && tall(j) >t5
                val = Xall(j-1);
                break;
            end
        end
        S5 = [S5, val];
    end
    % 计算衍生参数
    bs=nu/gamma; 
    bf=1/(1/lam1+1/lam2);
    % 并行临时数组赋值
    par_lam1(zushu)=lam1;
    par_lam2(zushu)=lam2;
    par_gamma(zushu)=gamma;
    par_nu(zushu)=nu;
    par_bf(zushu)=bf;
    par_bs(zushu)=bs;
    par_S5(zushu,:)=S5;
end
% 关闭并行池
delete(gcp);

% 回填结果到单元格数组
parameter(2:GN+1,2:5) = num2cell([par_lam1,par_lam2,par_gamma,par_nu]);
parameter(2:GN+1,7:8) = num2cell([par_bf,par_bs]);
expression(2:GN+1,2:N1+1) = num2cell(par_S5);

% 导出csv
writetable(cell2table(parameter), 'threestate_parameter500.csv');
writetable(cell2table(expression), 'threestate_expression500.csv');