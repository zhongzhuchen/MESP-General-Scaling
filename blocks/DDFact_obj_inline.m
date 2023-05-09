% obtain problem data
n=obj.size;
F=obj.F;
Fsquare = obj.Fsquare;
A_data=obj.A;
b_data=obj.b;
info = struct;
[fval,dx,~] = DDFact_obj_Knitro(x,s,F,Fsquare,Gamma);
% convert the sign back
dx=-dx;
fval=-fval;
%% calculate dual bound and dual solutions
% % cache for mixing
% info.cache1=-s;
% info.cache2=sum(dx)-s;
[m,~] = size(obj.A);

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
info.dualgap=dualgap-s+sum(x.*log(Gamma));
info.fval=fval;
info.dualbound=fval+info.dualgap;

% cache for mixing
info.cache = -s;
