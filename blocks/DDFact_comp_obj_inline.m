%% obtain problem data
info=struct;
n=obj.size;
A_data=obj.A;
b_data=obj.b;
[m,~] = size(obj.A);
F_comp=obj.F_comp;
Fsquare_comp = obj.Fsquare_comp;
ldetC = obj.ldetC;

[fval,dx,~] = DDFact_comp_obj_Knitro(x,s,F_comp,Fsquare_comp,ldetC,Gamma);
% convert the sign back
dx=-dx;
fval=-fval;
dx1 = dx+log(Gamma);
%% calculate dual solution
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
info.dualgap=dualgap-(n-s)-sum(dx)+sum((ones(n,1)-x).*log(Gamma));

%% transform objecvtive value and solution back
info.fval=fval;
info.dualbound=fval+info.dualgap;
info.cache=sum(dx1)-(n-s);
