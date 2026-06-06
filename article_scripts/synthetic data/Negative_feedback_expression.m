function SSA_negative_feedback_expression_matrix_parallel
clc
clear all
close all

% ===================== 核心配置 =====================
total_group = 10000;       % 基因
N = 4000;                % 细胞
TT_max = 200;           
H_time_step = 200;      
eps_steady = 0.01;      
% 输出文件
param_file = 'model_parameters4000.csv';       
expr_file = 'expression_matrix4000.csv';       
% ====================================================

% 初始化存储矩阵
param_matrix = zeros(total_group, 4);
expression_matrix = zeros(total_group, N);

% ===================== 并行循环核心修改 =====================
% parfor：自动调用Linux集群多核心，每个基因独立并行计算
parfor zushu = 1:total_group
    % 【删除了原代码的clearvars，并行下禁止使用！】
    fprintf('正在计算第 %d/%d 个基因...\n', zushu, total_group);

    % ========== 1. 随机参数采样（负反馈模型） ==========
    kon = 0.1 + 4.9*rand;
    koff = 0.1 + 2.9*rand;
    kb = 5 + 45*rand;
    nu = 0.05 + 0.95*rand;
    kd = 1;

    % ========== 2. 主方程求解稳态时间 ==========
    TT = TT_max; H = H_time_step; eps = eps_steady;
    TP = linspace(TT/H, TT, H);
    M = ceil(3*kb+1); dt = 0.001; S = TT/dt + 1;
    p1 = zeros(S,M+1); p2 = zeros(S,M+1);
    p1(1,1) = 1;

    for ta = 2:S
        p1(ta,1) = p1(ta-1,1) + dt*(kd*p1(ta-1,2) + koff*p2(ta-1,1) - kon*p1(ta-1,1));
        p2(ta,1) = p2(ta-1,1) + dt*(kd*p2(ta-1,2) + kon*p1(ta-1,1) - (kb+koff)*p2(ta-1,1));
        for i = 2:M
            p1(ta,i) = p1(ta-1,i) + dt*((koff + (i-1)*nu)*p2(ta-1,i) - ((i-1)*kd + kon)*p1(ta-1,i) + i*kd*p1(ta-1,i+1));
            p2(ta,i) = p2(ta-1,i) + dt*(kon*p1(ta-1,i) - (kb+(i-1)*(kd+nu)+koff)*p2(ta-1,i) + i*kd*p2(ta-1,i+1) + kb*p2(ta-1,i-1));
        end
    end
    pM = p1 + p2;

    meandata = zeros(1,H);
    seconddata = zeros(1,H);
    for j=1:H
        meandata(j) = sum((0:M) .* pM(round(TP(j)/dt)+1, 1:M+1));
        seconddata(j) = sum((0:M).^2 .* pM(round(TP(j)/dt)+1, 1:M+1));
    end
    kk = H; for j=1:H-1, if abs(meandata(H)-meandata(H-j))/meandata(H)<=eps, kk=kk-1; else, break; end; end
    kk1= H; for j=1:H-1, if abs(seconddata(H)-seconddata(H-j))/seconddata(H)<eps, kk1=kk1-1; else, break; end; end
    tend = max(kk*(TT/H), kk1*(TT/H));

    % ========== 3. SSA随机模拟 ==========
    T = tend;
    expr = zeros(1, N);

    for sim = 1:N
        X = 0; t = 0; s = 1;
        while t <= T
            h4 = X;
            a1 = kb;
            a2 = koff + nu*X;
            a3 = kon;
            a4 = kd*X;

            if s == 0
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
            else
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
        expr(sim) = X;
    end

    % 并行赋值（parfor标准用法，安全无冲突）
    param_matrix(zushu, :) = [kon, koff, kb, nu];
    expression_matrix(zushu, :) = expr;
end
% ========================================================

% ========== 写入CSV文件（串行写入，避免文件冲突） ==========
fprintf('\n正在写入文件...\n');

fid = fopen(param_file,'w');
fprintf(fid,'gene_id,kon,koff,kb,nu\n');
for i=1:total_group
    fprintf(fid,'gene%d,%.6f,%.6f,%.6f,%.6f\n',i,param_matrix(i,1),param_matrix(i,2),param_matrix(i,3),param_matrix(i,4));
end
fclose(fid);

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

fprintf('\n✅ 并行计算完成！\n');
end