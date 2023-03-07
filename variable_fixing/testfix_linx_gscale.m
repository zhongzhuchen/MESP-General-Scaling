clear all;
load('data124.mat');
n = length(C);
%% open parameter
maxfixon = 0;
multifix_threshold0 = 0.05;
multifix_threshold1 = 0.95;
A = double.empty(0,n);
b = double.empty(0,1);
% load('linconstr124.mat')
MESPEX1 = MESP(C, A, b);
exceloutput=[];
folder = 'variable_fixing/results_constr';
ubname = 'linx';

baseFileName = strcat(ubname, '_data',int2str(n),'.xlsx');
fullFileNameexcel = fullfile(folder, baseFileName);
if exist(fullFileNameexcel, 'file')==2
    delete(fullFileNameexcel);
end

if maxfixon == 0
    info1.fix0num = nan;
    info1.fix1num = nan;
    info1.fixnum = nan;
end

% collect the variables fixed for each s
info_fixconflict = struct;
info_heavyfix = struct;
info_maxfix = struct;
info_maxfixlin = struct;
info_fix = struct;

if ubname == 'linx'
    for s=11:110
        sprintf("s = %d", s)
        suboutput = [n, s];
        LB = MESPEX1.obtain_lb(s);

        % g-scaling
        [optgamma,info]= MESPEX1.Newton_Linx_gamma(s, 300);
        [optGamma,info]= MESPEX1.BFGS_Linx_Gamma(s, sqrt(optgamma)*ones(n,1), 300);
        [fval,x,info] = MESPEX1.Knitro_Linx(s/n*ones(n,1),s,optGamma);
        [fval,x,info0] = MESPEX1.Knitro_Linx_heavy(s/n*ones(n,1),s,optGamma);
        if maxfixon == 1
            [fix0vec, fix1vec, info1] = maxfix(ubname, 1:n, C, s, A, b, s/n*ones(n,1), optGamma, LB);
        end
        [fix0vec, fix1vec, info2] = maxfixlin(ubname, 1:n, C, s, A, b, s/n*ones(n,1), optGamma, LB);
        [fix0vec, fix1vec, info3] = singlefix(ubname, 1:n, C, s, A, b, s/n*ones(n,1), optGamma, LB);
        
        % multifix: rule 1
        I1 = [];
        J1 = [];
        for i=1:n
            if info.x(i) < multifix_threshold0
                I1(end+1) = i;
            elseif info.x(i) > multifix_threshold1
                J1(end+1) = i;
            end
        end
        [fix0vec, fix1vec, info4] = multifix(ubname, C, s, A, b, I1, J1, s/n*ones(n,1), optGamma, LB);

        % multifix: rule 2
        I2 = [];
        J2 = [];
        for i=1:n
            if info.fval-LB-info.dual_upsilon(i) < multifix_threshold0
                I2(end+1) = i;
            elseif info.fval-LB-info.dual_nu(i) < multifix_threshold0
                J2(end+1) = i;
            end
        end
        [fix0vec, fix1vec, info5] = multifix(ubname, C, s, A, b, I2, J2, s/n*ones(n,1), optGamma, LB);

        suboutput = [suboutput, info.solved, ...
            info.fix0num, info0.fix0num, info1.fix0num, info2.fix0num, info3.fix0num, info4.fix0num, info5.fix0num,...
            info.fix1num, info0.fix1num, info1.fix1num, info2.fix1num, info3.fix1num, info4.fix1num, info5.fix1num,...
            info.fixnum, info0.fixnum, info1.fixnum, info2.fixnum, info3.fixnum, info4.fixnum, info5.fixnum];

        exceloutput(end+1,:) = suboutput;
    end
    title=["n", "s", 'solved',...
        'g-scale CG 0', 'g-scale heavyfix 0', 'g-scale maxfix 0', 'g-scale maxfixlin 0', 'g-scale singlefix 0', 'g-scale multifix 0', 'g-scale multifix rule 2 0',...
        'g-scale CG 1', 'g-scale heavyfix 1', 'g-scale maxfix 1', 'g-scale maxfixlin 1', 'g-scale singlefix 1', 'g-scale multifix 1', 'g-scale multifix rule 2 1',...
        'g-scale CG all', 'g-scale heavyfix all', 'g-scale maxfix all', 'g-scale maxfixlin all', 'g-scale singlefix all','g-scale multifix all', 'g-scale multifix rule 2 all'];
end
xlswrite(fullFileNameexcel ,title,1,'A1');
xlswrite(fullFileNameexcel ,exceloutput,1,'A2');