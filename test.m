%% test new heur method
load('data63.mat');
n=length(C);
m=5;
A = double.empty(0,n);
b = double.empty(0,1);
% load('random_lincon.mat');
MESPInstance1 = MESP(C,A,b);
% LB = [];
% Heurval =[];
% for s=10% 10:(n-1)
%     s
%     LB(end+1)=a.obtain_lb(s);
%     [xind, heurval] = heur(C,s,A,b);
%     A1 = randi([1,10],m,n);
%     x=zeros(n,1);
%     x(xind)=1;
%     b1 = A1*x-ones(m,1);
%     MESPInstance1 = MESP(C,A1,b1);
%     Heurval(end+1)=MESPInstance1.obtain_lb(s);
% end
% compare = [LB', Heurval'];

%%test feasible
% for s=10:53
%     [lb,info] = Prob1.obtain_lb(s);
% end

% allx=zeros(n,length(2:58));
% for s=2:58
%     x0=s/n*ones(n,1);
% %     [gamma,~]=Prob.BFGS_Linx_gamma(s);
%     [fval,x,info]=Prob.Knitro_DDFact(x0,s,ones(n,1));
%     allx(:,s-1)=x;
% end

% clear all;
% load('data63.mat');
% systermatic_experiments_scalar2;
% systermatic_experiments_scalar3;
% clear all;
% load('data90.mat');
% systermatic_experiments_scalar;
% clear all;
% load('data124.mat');
% systermatic_experiments_scalar;

% for s=5:5:50
%     [fval,x,info] = Prob.Knitro_DDFact_comp(s/n*ones(n,1),s,ones(n,1));time=info.time;
%     [fval,dx,info] = Prob.DDFact_comp_obj(x,s,ones(n,1));
%     info.dualtime/time
% end

s = 5;
Gamma = ones(n,1);
% [fval,x,info] = SDPT3_BQP_light(diag(s/n*ones(n,1)),C,s,A,b,Gamma);
% [fval,x,info] = MESPInstance.SDPT3_BQP(diag(s/n*ones(n,1)), s, Gamma);
% [optGamma,info] = MESPInstance.BFGS_BQP_Gamma(s,Gamma);
