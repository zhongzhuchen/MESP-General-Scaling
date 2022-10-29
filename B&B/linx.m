function [results,xval,delta_zero,delta_one]=linx(C,s,control)
%
% evaluate linx bound for MESP. Return bound, primal and dual variables,
% updated gamma value and stats. Can also evaluate dual bound if desired.

n=size(C,1);
%
% nscale=control(1);
% nscale0=control(2);
% nrounds=control(3);   % not used by linx 
% Maxadd=control(4);    % not used by linx
% printsol=control(5);
% gamma=control(6);
% power=control(7); 
% maxfactor=control(8);

nscale = control.nscale;
nscale0 = control.nscale0;
nrounds = control.nrounds; % not used by linx 
Maxadd = control.Maxadd; % not used by linx 
printsol = control.printsol;
Gamma = (control.Gamma)';
power = control.power;
maxfactor = control.maxfactor;

%    
if Gamma==0
%     csort=sort(diag(C),'descend');
%     gamma=1/csort(s); % initial scale factor if none provided
%     gamma=gamma^power;
    A = double.empty(0,n);
    b = double.empty(0,1); 
    MESPInstance = MESP(C,A,b);
    [optgamma,info_Linxg1]= MESPInstance.Newton_Linx_gamma(s, 300);
    [Gamma,info_Linxg2]= MESPInstance.BFGS_Linx_Gamma(s, sqrt(optgamma)*ones(n,1), 300);
    nscale=nscale0;
end
%
x=sdpvar(n,1);
%
time=0; % total time spent in solvers
%
for kscale=1:nscale    
    %
    if kscale>1; Gamma=newGamma; end % apply adjustment from previous iteration
    %
    constraints=[x >= zeros(n,1), x <= ones(n,1), sum(x)==s];
    scaleC = diag(Gamma)*C;
    A = scaleC*diag(x)*scaleC;
    obj1=0.5*(logdet(.5*(A+A')+eye(n)-diag(x)))-sum(x.*log(Gamma)); % force symmetry
    options=sdpsettings('solver','sdpt3','sdpt3.maxit',75,'sdpt3.gaptol',1E-7,'sdpt3.inftol',1E-7,'verbose',0,'cachesolvers',1);
    solvetime=tic;
    diagnostics=optimize(constraints,-obj1,options);    
    time=time+toc(solvetime); % record time in solver
    code=diagnostics.problem;
    %  
    xval=value(x);
    AUX=scaleC*diag(xval)*scaleC';
    F=AUX+eye(n)-diag(xval);
    F=(F+F')/2; % force symmetry
%     bound=.5*(logdet(F)-s*log(gamma)); % Note factor .5: bound is for 2*ldet C_FF
    bound = double(obj1);
    %
    [U,D]=eig(F);
    lam=diag(D);
    Finv=U*diag(1./lam)*U';
    scaleval=norm(diag(Finv*AUX)-xval); % want scaleval = 0
    if printsol; fprintf('\n code=%g, bound=%g', code,bound); end
    if printsol; fprintf('\n gamma = %g, n-s= %g, scaleval = %g',Gamma, n-s, scaleval); end
    newGamma=Gamma;
    scalerror=abs(scaleval);
    if scalerror <.1; break; end % don't even try to adjust once close enough
%     deriv=diag(Finv*C*diag(xval)*C*Finv)'*(xval-ones(n,1));
%     newgamma=gamma+.8*((n-s)-scaleval)/deriv; %damp step
%     newgamma=min(newgamma,maxfactor*gamma); % enforce max factor change on one step
%     newgamma=max(newgamma,gamma/maxfactor);
    %
    linx_Gamma_update;
    if scalerror <.5; break; end % good enough for current bound
end
%
delta_one=-.5*dual(constraints(1));  % change sign and correct for factor 2 in objective
delta_zero=-.5*dual(constraints(2)); %
%
% evaluate_dual=0; % to evaluate dual bound
% % fixme: modify evaluate dual
% if evaluate_dual;
%     %
%     S=sdpvar(n,n);
%     v=sdpvar(n,1);
%     u=sdpvar(1,1);
%     %
%     dual_constraints=[ S >= 0, v >= 0, gamma*diag(C*S*C) - diag(S) <= u*ones(n,1)+v];
%     dual_constraints=[dual_constraints, trace(S)+ s*u + ones(1,n)*v == n];
%     diagnostics=optimize(dual_constraints,-logdet(S),options);    
%     code=diagnostics.problem;
%     %
%     Sval=value(S);
%     dualbound=.5*(-logdet(Sval)-s*log(gamma));
%     vval=value(v);
%     uval=value(u);
%     wval=uval*ones(n,1)+vval+diag(Sval)-gamma*diag(C*Sval*C);
%     if printsol; fprintf('\n code=%g, dual bound=%g', code,dualbound); end
%     % [xval wval vval delta_one delta_zero]
% end
%
results = struct;
results.bound=bound;
results.code=code;
results.time=time;
results.newGamma=newGamma';
results.kscale=kscale;
results.scalerror=scalerror;

