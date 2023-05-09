function [fval,x,info] = SDPT3_BQP_comp_light(Y0,Cinv,s,A,b,ldetC,Gamma)
%% obtain class properties
[m,~] = size(A);
scaleCinv=diag(Gamma)*Cinv*diag(Gamma);
n = length(Cinv);
info = struct;
%% calling SDPT3
Y=sdpvar(n,n);
assign(Y,Y0);
x=ones(n,1)-diag(Y);  
Y1=[Y, ones(n,1)-x; (ones(n,1)-x)',1];  
%constraints=[ Y >=0, sum(x)==s, sum(X)==s*x'];
% or use the following if numerical trouble:
% constraints=[ Y >=0, 1.0000001*s >= sum(x)>=0.9999999*s,  1.0000001*s*x' >= sum(X)>=0.9999999*s*x', A*x-b <= 0];
constraints=[ Y1 >=0, sum(x)==s,  sum(Y)==(n-s)*(ones(n,1)-x)', A*x-b <= 0];
options=sdpsettings('solver','sdpt3','sdpt3.maxit',100,'sdpt3.gaptol',1E-10,...
    'sdpt3.inftol',1E-10,'sdpt3.steptol',1E-10,'verbose',0);
obj1= logdet(scaleCinv.*Y+diag(x))-2*sum((ones(n,1)-x).*log(Gamma));
TStart=tic;
tStart=cputime;
diagnostics=optimize(constraints,-obj1,options);
time=toc(TStart);
tEnd=cputime-tStart;

Yval=value(Y);
x=ones(n,1)-diag(Yval);
% ========================================
fval=double(obj1)+ldetC;
info.exitflag=diagnostics.problem;
info.fval=fval;
info.Y=Yval;
info.x=x;
info.time=time;
info.cputime = tEnd;
info.info = diagnostics.info;
end
