function [fval,x,info] = SDPT3_BQP_light(X0,C,s,A,b,Gamma)
%% obtain class properties
[m,~] = size(A);
scaleC=diag(Gamma)*C*diag(Gamma);
n = length(C);
info = struct;
%% calling SDPT3
X=sdpvar(n,n);
assign(X,X0);
x=diag(X);  
Y=[X, x; x',1];  
%constraints=[ Y >=0, sum(x)==s, sum(X)==s*x'];
% or use the following if numerical trouble:
% constraints=[ Y >=0, 1.000001*s >= sum(x)>=0.999999*s,  1.000001*s*x' >= sum(X)>=0.999999*s*x', A*x-b <= 0];
constraints=[ Y >=0, sum(x)==s,  sum(X)==s*x', A*x-b <= 0];
options=sdpsettings('solver','sdpt3','sdpt3.maxit',100,'sdpt3.gaptol',1E-10,...
    'sdpt3.inftol',1E-10,'sdpt3.steptol',1E-10,'verbose',0);
obj1= logdet(scaleC.*X+eye(n)-diag(x))-2*sum(x.*log(Gamma));
TStart=tic;
tStart=cputime;
diagnostics=optimize(constraints,-obj1,options);
time=toc(TStart);
tEnd=cputime-tStart;

Xval=value(X);
x=diag(Xval);
% ========================================
fval=double(obj1);
info.exitflag=diagnostics.problem;
info.fval=fval;
info.X=Xval;
info.x=x;
info.time=time;
info.cputime = tEnd;
info.info = diagnostics.info;
end
