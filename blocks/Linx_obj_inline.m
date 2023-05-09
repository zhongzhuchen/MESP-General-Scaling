%% obtain problem data
C=obj.C;
n=obj.size;
A_data=obj.A;
b_data=obj.b;
[m,~] = size(obj.A);
info = struct;
%% calculate the objective value and gradient
[fval,dx,ininfo] = Linx_obj_Knitro(x,C,Gamma);
% convert the sign back
dx=-dx;
fval=-fval;
dx2=ininfo.dx2;
%% calculate dual bound and dual solutions
f=[zeros(n,1);ones(n,1);b_data;s];
Aeq=[-eye(n),eye(n),A_data',ones(n,1)];
beq=dx;
lb=[zeros(2*n+m,1);-inf];
ub=Inf(2*n+m+1,1);
x0=[];
options = knitro_options('algorithm',3,...  % active-set/simplex algorithm
                         'outlev',0);       % iteration display
tstart=tic;
[xlp, dualgap, exitflag, ~] = knitro_lp (f, [], [], Aeq, beq, lb, ub, x0, [], options);
time=toc(tstart);
info.dualtime=time;

info.dual_upsilon = xlp(1:n);
info.dual_nu = xlp((n+1):2*n);
info.dual_pi =  xlp((2*n+1):(2*n+m));
info.dual_tau = xlp(end);
info.dualgap=dualgap+0.5*sum(dx2)-n/2+sum(x.*log(Gamma));
info.fval=fval;
info.dualbound=fval+info.dualgap;
% cache for mixing
info.cache = 0.5*sum(dx2)-n/2;


