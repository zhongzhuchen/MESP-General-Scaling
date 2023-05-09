n=3;
s=2;

tau=sdpvar(1,1);
pi=sdpvar(1,1);
upsilon=sdpvar(n,1);
nu=sdpvar(n,1);

epsilon = 0.25;
c=[-1;-2-epsilon;-3];
b=1.5;
time=0; % total time spent in solvers


constraints=[tau*ones(n,1)+pi*[ones(2,1);0]+nu-upsilon==c,pi >= 0,nu>= zeros(n,1),upsilon >= zeros(n,1)];
%Either use this:
% obj1= s*tau+1.5*pi+sum(nu)+ 4 - nu(1);
%
%or this:
obj1=0;
constraints=[constraints,s*tau+b*pi+sum(nu)+ 4 - nu(1)<=1e-6,s*tau+b*pi+sum(nu)+ 4 - upsilon(2)<=1e-6];
%
options=sdpsettings('solver','sdpt3','sdpt3.maxit',75,'sdpt3.gaptol',1E-7,'sdpt3.inftol',1E-7,'verbose',1,'cachesolvers',1);
solvetime=tic;
diagnostics=optimize(constraints,obj1,options);    
time=time+toc(solvetime) % record time in solver
code=diagnostics.problem
%  
tauval=value(tau)
pival=value(pi)
upsilonval=value(upsilon)
nuval=value(nu)
objval=value(obj1)

