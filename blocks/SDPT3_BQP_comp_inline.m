%% obtain problem data
C=obj.C;
Cinv=obj.C_comp;
ldetC=obj.ldetC;
n=obj.size;
A_data=obj.A;
b_data=obj.b;
A_comp_data = -A_data;
b_comp_data = b_data-A_data*ones(n,1);
[m,~] = size(obj.A);
info = struct;

%% calling SDPT3
[fval,x,ininfo] = SDPT3_BQP_comp_light(Y0,Cinv,s,A_data,b_data,ldetC,Gamma);
%% assign values to info
% record important information
info.exitflag=ininfo.exitflag;
info.Y=ininfo.Y;
info.x=x;

info.fval=fval;
info.time=ininfo.time;
info.cputime=ininfo.cputime;

%% calculate dualgap and bound
% % % Y=ininfo.Y;
x=ininfo.x;
Y=ininfo.Y;
Ytilde = [Y, ones(n,1)-x; (ones(n,1)-x)', 1];
scaleCinv = diag(Gamma)*Cinv*diag(Gamma);
Ctilde = [scaleCinv-eye(n), zeros(n,1); zeros(1,n), 0];
Beta = Ctilde.*Ytilde+eye(n+1);
Theta = inv(Beta);

% we formulate the semidefinite programming problem for calculating dualgap
% as follows
pi = sdpvar(m,1);
tau = sdpvar(2*n+2,1);
constformula = 0;
h1 = b_comp_data;
h2 = zeros(2*n+2,1);
for i = 1:m
    constformula = constformula + pi(i)*[zeros(n), 0.5*A_comp_data(i,:)'; 0.5*A_comp_data(i,:),0];
end
I = eye(n);
for i=1:n
    constformula = constformula + tau(i)*[I(:,i)*I(:,i)', -0.5*I(:,i); -0.5*I(:,i)', 0];
end
for i=1:n
    Ki = 0.5*(diag(I(:,i))*ones(n)+(diag(I(:,i))*ones(n))');
    constformula = constformula + tau(i+n)*[Ki, -0.5*(n-s)*I(:,i); -0.5*(n-s)*I(:,i)', 0];
end
constformula = constformula + tau(2*n+1)*[zeros(n), 0.5*ones(n,1); 0.5*ones(1,n), 0];
h2(2*n+1)=(n-s);
constformula = constformula + tau(2*n+2)*[zeros(n), zeros(n,1); zeros(1,n), 1];
h2(2*n+2)=1;

constformula = constformula + [zeros(n), log(Gamma); log(Gamma)', 0];

constraints = [constformula >= Theta.*Ctilde; pi>=0];
options=sdpsettings('solver','sdpt3','sdpt3.maxit',100,'sdpt3.gaptol',1E-10,...
    'sdpt3.inftol',1E-10,'sdpt3.steptol',1E-10,'verbose',0);
obj1 = pi'*h1+(n-s)*tau(2*n+1)+tau(2*n+2);
TStart=tic;
tStart=cputime;
diagnostics=optimize(constraints,obj1,options);
time=toc(TStart);
tEnd=cputime-tStart;
dualgap = double(obj1)+trace(Theta)-(n+1)+2*sum((ones(n,1)-x).*log(Gamma));
info.continuous_dualgap = dualgap;
if dualgap > 0
    info.dualbound = fval+dualgap;
else
    info.dualbound = fval;
end

%% calculate upsilon/nu in MESP book Theorem 3.6.12
upsilon = zeros(n,1);
nu = zeros(n,1);
for i=1:n
    if scaleCinv(i,i) <= 0
        upsilon(i) = inf;
        nu(i)= -inf;
    else
        upsilon(i) = (scaleCinv(i,i)+1)/Beta(i,i) - log((scaleCinv(i,i)+1)/Beta(i,i)) -1;
        nu(i) = 1/Beta(i,i) - log(1/Beta(i,i)) -1;
    end
end
% exchange
info.dual_upsilon = nu;
info.dual_nu = upsilon;

%% fixing variables by the fixing logic in the factorization paper by Chen etal.
info.fixnum=0;
info.fixnum_to0=0;
info.fixto0list=[];
info.fixnum_to1=0;
info.fixto1list=[];

% if the integrality gap between the upper bound and the lower bound is
% large, the fixing power might be hidden due the weak lower bound
info.integrality_gap=info.dualbound-obj.obtain_lb(s);
if info.integrality_gap>1e-6
    info.solved=0;
else
    info.solved=1;
end
for i=1:n
    if info.integrality_gap<info.dual_upsilon(i)-1e-10 % fix to zero
        info.fixnum=info.fixnum+1;
        info.fixnum_to0=info.fixnum_to0+1;
        info.fixto0list(end+1)=i;
    elseif info.integrality_gap<info.dual_nu(i)-1e-10 % fix to one
        info.fixnum=info.fixnum+1;
        info.fixnum_to1=info.fixnum_to1+1;
        info.fixto1list(end+1)=i;
    end
end

