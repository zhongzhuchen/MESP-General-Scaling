function [fval,x,info] = Knitro_DDFact_comp_light_check(x0,C,s,F_comp,Fsquare_comp,ldetC,A,b,Gamma,fix)
%% obtain class properties
A_data=A;
b_data=b;
info=struct;
n =length(x0);
%% calling knitro
obj_fn =  @(x) DDFact_comp_obj_Knitro(x,s,F_comp,Fsquare_comp,ldetC,Gamma);
lb=zeros(n,1);
ub=ones(n,1);
vfix0 = fix.fixto0list;
vfix1 = fix.fixto1list;
IM = eye(n);
Aeq=[ones(1,n); IM(vfix1, :); IM(vfix0, :)];
beq=[s;ones(length(vfix1),1);zeros(length(vfix0),1)];
A=A_data;
b=b_data;

if sum(abs(Aeq*x0-beq))>n*1e-10
    error('The initial point x0 is not feasible.')
end

options = knitro_options('algorithm', 0, 'convex', 1, 'derivcheck', 0, 'outlev', 0 , 'gradopt', 1, ...
                         'hessopt', 2, 'maxit', 1000, 'xtol', 1e-15, ...
                         'feastol', 1e-10, 'opttol', 1e-10, 'bar_feasible',1,...
                         'bar_maxcrossit', 10);
TStart=tic;
tStart=cputime;
[x,knitro_fval,exitflag,output,lambda,~] = knitro_nlp(obj_fn,x0,A,b,Aeq,beq,lb,ub,[],[],options);  
time=toc(TStart);
tEnd=cputime-tStart;

fval=-knitro_fval;
info.exitflag=exitflag;
info.output=output;
info.lambda=lambda;
info.fval=fval;
info.time=time;
info.cputime = tEnd;
end

