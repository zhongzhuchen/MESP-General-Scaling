n=3;
s=2;
c=[-1;-5/2;-3];

% Solving the primal:
x=sdpvar(n,1);
time=0; % total time spent in solvers
    constraints=[x >= zeros(n,1), x <= ones(n,1), sum(x)==s,x(1)+x(2)<=1.5];
    obj1= -x(1)-2.5*x(2)-3*x(3);
    options=sdpsettings('solver','sdpt3','sdpt3.maxit',75,'sdpt3.gaptol',1E-7,'sdpt3.inftol',1E-7,'verbose',0,'cachesolvers',1);
    solvetime=tic;
    diagnostics=optimize(constraints,-obj1,options);    
    time=time+toc(solvetime) % record time in solver
    code=diagnostics.problem
    xval=value(x)
    objval=value(obj1)
    upsilonval1=dual(constraints(1))  % change sign and correct for factor 2 in objective
    nuval1=dual(constraints(2)) 

%Solving the dual:
tau=sdpvar(1,1);
pi=sdpvar(1,1);
upsilon=sdpvar(n,1);
nu=sdpvar(n,1);
timedual=0; % total time spent in solvers
    constraints=[tau*ones(n,1)+pi*[ones(2,1);0]-upsilon+nu==c,pi >= zeros(n,1),upsilon >= zeros(n,1),nu>= zeros(n,1) ];
    obj1= s*tau+1.5*pi+sum(nu);
    options=sdpsettings('solver','sdpt3','sdpt3.maxit',75,'sdpt3.gaptol',1E-7,'sdpt3.inftol',1E-7,'verbose',0,'cachesolvers',1);
    solvetime=tic;
    diagnostics=optimize(constraints,obj1,options);    
    timedual=timedual+toc(solvetime) % record time in solver
    code=diagnostics.problem
    tauval=value(tau)
    pival=value(pi)
    upsilonval=value(upsilon)
    nuval=value(nu)
    objval=value(obj1)
   
