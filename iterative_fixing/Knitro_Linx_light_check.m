function [fval,x,info] = Knitro_Linx_light_check(x0,C,s,A,b,Gamma,fix)
%% obtain class properties
A_data=A;
b_data=b;
[m,~] = size(A);
scaleC=diag(Gamma)*C;
n = length(x0);
info = struct;
%% calling knitro
obj_fn =  @(x) Linx_obj_Knitro(x,C,Gamma);
lb=zeros(n,1);
ub=ones(n,1);
vfix0 = fix.fixto0list;
vfix1 = fix.fixto1list;
IM = eye(n);
Aeq=[ones(1,n); IM(vfix1, :); IM(vfix0, :)];
beq=[s;ones(length(vfix1),1);zeros(length(vfix0),1)];
A=A_data;
b=b_data;

% if sum(abs(Aeq*x0-beq))>n*1e-8
%     error('The initial point x0 is not feasible.')
% end

TStart=tic;
tStart=cputime;

extendedFeatures.HessFcn = @(x,lambda) Linx_hessfun(x,lambda,scaleC);
options = knitro_options('algorithm', 3, 'convex', 1, 'derivcheck', 0, 'outlev', 0 , 'gradopt', 1, ...
                         'hessopt', 1, 'maxit', 1000, 'xtol', 1e-15, 'derivcheck_tol',1e-5,...
                         'feastol', 1e-10, 'opttol', 1e-10, 'bar_feasible',1,...
                         'bar_maxcrossit', 10);
[x,knitro_fval,exitflag,output,lambda,~] = knitro_nlp(obj_fn,x0,A,b,Aeq,beq,lb,ub,[],extendedFeatures,options);
time=toc(TStart);
tEnd=cputime-tStart;
% ========================================
fval=-knitro_fval;
info.exitflag=exitflag;
info.output=output;
info.lambda=lambda;
info.fval=fval;
info.time=time;
info.cputime = tEnd;
end
