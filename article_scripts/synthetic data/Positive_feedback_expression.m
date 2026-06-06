function SSA_positive_feedback_expression_matrix_parallel
clc
clear all
close all

% ===================== 核心配置（与负反馈代码格式一致） =====================
total_group = 10000;         % 基因总数（原代码500个）
N = 500;                   % 细胞数
TT_max = 200;              % 主方程最大计算时间
H_time_step = 200;         % 主方程时间分段数
eps_steady = 0.01;         % 稳态判断阈值
% 输出CSV文件（格式与负反馈完全相同）
param_file = 'positive_model_parameters500.csv';       
expr_file = 'positive_expression_matrix500.csv';       
% ============================================================================

% 初始化存储矩阵（与负反馈格式一致）
param_matrix = zeros(total_group, 4);
expression_matrix = zeros(total_group, N);

% ===================== 并行循环（与负反馈完全一致） =====================
parfor zushu = 1:total_group
    fprintf('正在计算正反馈模型 第 %d/%d 个基因...\n', zushu, total_group);

    % ========== 1. 正反馈模型随机参数采样（原代码参数范围，保留不变） ==========
    kon  = 0.1 + 2.9*rand;   % 基因激活速率
    koff = 0.1 + 4.9*rand;   % 基因失活速率
    kb   = 5 + 45*rand;      % 转录速率
    mu   = 0.05 + 1.95*rand; % 正反馈调控强度（核心参数）
    kd   = 1;                % 降解速率（固定）

    % ========== 2. 化学主方程求解稳态时间（正反馈公式，保留不变） ==========
    TT = TT_max; H = H_time_step; eps = eps_steady;
    TP = linspace(TT/H, TT, H);
    M = ceil(3*kb+1);
    dt = 0.001;
    S = TT/dt + 1;
    
    % 概率矩阵初始化：p1=沉默态，p2=激活态
    p1 = zeros(S,M+1); p2 = zeros(S,M+1);
    p1(1,1) = 1; p2(1,1) = 0;

    % 主方程迭代（正反馈核心方程，未修改）
    for ta = 2:S 
        p1(ta,1) = p1(ta-1,1) + dt*(kd*p1(ta-1,2) + koff*p2(ta-1,1) - kon*p1(ta-1,1));
        p2(ta,1) = p2(ta-1,1) + dt*(kd*p2(ta-1,2) + kon*p1(ta-1,1) - (kb+koff)*p2(ta-1,1));
        for i = 2:M  
            p1(ta,i) = p1(ta-1,i) + dt*(koff*p2(ta-1,i) - ((i-1)*(kd+mu)+kon)*p1(ta-1,i) + i*kd*p1(ta-1,i+1));
            p2(ta,i) = p2(ta-1,i) + dt*((kon+(i-1)*mu)*p1(ta-1,i) - (kb+(i-1)*kd+koff)*p2(ta-1,i) + i*kd*p2(ta-1,i+1) + kb*p2(ta-1,i-1));
        end
    end
    pM = p1 + p2; % 总概率分布

    % 计算均值、二阶矩，判断稳态
    meandata = zeros(1,H);
    seconddata = zeros(1,H);
    for j=1:H
        meandata(j) = sum((0:M) .* pM(round(TP(j)/dt)+1, 1:M+1));
        seconddata(j) = sum((0:M).^2 .* pM(round(TP(j)/dt)+1, 1:M+1));
    end

    % 稳态判断逻辑（与负反馈完全一致）
    kk = H; 
    for j=1:H-1
        if abs(meandata(H)-meandata(H-j))/meandata(H) <= eps
            kk = kk-1;
        else
            break;
        end
    end
    kk1 = H; 
    for j=1:H-1
        if abs(seconddata(H)-seconddata(H-j))/seconddata(H) < eps
            kk1 = kk1-1;
        else
            break;
        end
    end
    tend = max(kk*(TT/H), kk1*(TT/H)); % 最终稳态时间

    % ========== 3. SSA随机模拟（正反馈模型，与负反馈逻辑对齐） ==========
    T = tend;
    expr = zeros(1, N); % 存储当前基因100个细胞的表达量

    for sim = 1:N
        X = 0;    % 初始分子数
        t = 0;    % 初始时间
        s = 1;    % 初始状态：1=沉默，0=激活

        while t <= T
            % 正反馈反应速率（原代码核心逻辑不变）
            a1 = kb;          % 转录生成
            a2 = koff;        % 激活→沉默
            a3 = kon + mu*X;  % 正反馈核心：沉默→激活（分子越多，激活越快）
            a4 = kd*X;        % 分子降解

            if s == 0  % 基因处于激活态
                a0 = a1 + a2 + a4;
                tau = -log(rand)/a0;
                r = rand*a0;
                if r <= a1
                    X = X + 1;
                elseif r <= a1+a2
                    s = 1;
                else
                    X = X - 1;
                end
            else  % 基因处于沉默态
                a0 = a3 + a4;
                tau = -log(rand)/a0;
                r = rand*a0;
                if r <= a3
                    s = 0;
                else
                    X = X - 1;
                end
            end
            t = t + tau;
        end
        expr(sim) = X; % 记录稳态表达量
    end

    % 并行安全赋值（与负反馈完全一致）
    param_matrix(zushu, :) = [kon, koff, kb, mu];
    expression_matrix(zushu, :) = expr;
end
% ========================================================================

% ========== 串行写入CSV文件（格式与负反馈100%相同，无文件冲突） ==========
fprintf('\n正在写入正反馈模型结果文件...\n');

% 1. 写入参数文件：gene_id,kon,koff,kb,mu
fid = fopen(param_file,'w');
fprintf(fid,'gene_id,kon,koff,kb,mu\n');
for i=1:total_group
    fprintf(fid,'gene%d,%.6f,%.6f,%.6f,%.6f\n',i,param_matrix(i,1),param_matrix(i,2),param_matrix(i,3),param_matrix(i,4));
end
fclose(fid);

% 2. 写入表达矩阵：gene_id + 100个细胞的表达量
fid = fopen(expr_file,'w');
fprintf(fid,'gene_id');
for i=1:N, fprintf(fid,',cell%d',i); end
fprintf(fid,'\n');
for i=1:total_group
    fprintf(fid,'gene%d',i);
    for j=1:N, fprintf(fid,',%d',expression_matrix(i,j)); end
    fprintf(fid,'\n');
end
fclose(fid);

fprintf('\n✅ 正反馈模型并行计算完成！\n');
end