%% obtain class properties
A_data=obj.A;
b_data=obj.b;
[m,~] = size(obj.A);
info=struct;
C = obj.C;
n = obj.size;
t1=tic;

if n>200
    TOL= 10^(-6);
    Numiterations=20; 
else
    TOL= 10^(-10);
    Numiterations=200; 
end

%% calculate the better lower bound among C and Cinv if C is invertible
heurval = obj.obtain_lb(s);

X0=diag(s/n*ones(n,1));
% Initialize Gamma
Gamma=GammaInit;

difgap=1;
k=1;

c1=1e-4;
c2=0.95;
% timelimit = 350;

%% calculate the gradient of linx bound with respect to Gamma
[bound,x,ininfo] = SDPT3_BQP_light(X0,C,s,A_data,b_data,Gamma);

X=ininfo.X;
scaleC=diag(Gamma)*C*diag(Gamma);
AUX = scaleC.*X;
B = - diag(x) + eye(n);
%Compute F(gamma,x)
F= AUX+B;
F=(F+F')/2; % force symmetry
%Compute inv(F(gamma,x))
[U,D]=eig(F);
lam=diag(D);
Finv=U*diag(1./lam)*U';
%Compute the residual res
gap=bound-heurval;

grad=2*(diag(Finv*AUX)-x);
res= norm(grad);

allGamma=Gamma;
allres=res;
allbound=bound;

H=eye(n); % initialize the inverse Hessian approximation

%% we use the optimal solution of every last linx bound as the initial point 
% for solving the next linx bound, trick for accelarating optimization
nX=X;
%% loop
while(k<=Numiterations && gap > TOL && abs(res) > TOL && difgap > TOL && toc(t1)<= timelimit)
    sprintf('iteration: %d, res: %f',k,abs(res));
    if k>1
        difgap=abs(allbound(k)-allbound(k-1));
        if k>=2 
            s_gmuly=deltag'*deltay;
            if k==2 
                gam=s_gmuly/(deltag'*deltag);
                H=gam*H;
            end
            if(s_gmuly>=TOL)
                Hg=H*deltag;
                Hgy=Hg*deltay';
                gHg=deltag'*Hg;
                yg=s_gmuly;
                H=H+(yg+gHg)/(yg^2)*(deltay*deltay')-(Hgy+Hgy')/yg;
            end
        end
    end
    if rank(H)<n
        H=eye(n);
    end
    %Compute the search direction for Newton's method (dir=-res/Mat)
    dir=-H*grad;
    edir=exp(dir);
    if norm(edir)==Inf || norm(edir)==0 || isnan(norm(edir))
        if norm(dir)==Inf || isnan(norm(dir))
            dir=-grad;
        else
            dir=dir/norm(dir)*10;
        end   
    end
    %check if alfa=1 satisfies the Strong Wolfe Conditions
    alfa=1;
    nGamma=Gamma.*exp(alfa*dir);
    [nbound,nx,ininfo] = SDPT3_BQP_light(nX,C,s,A_data,b_data,nGamma);
    nX=ininfo.X;
    scaleC=diag(nGamma)*C*diag(nGamma);
    nAUX = scaleC.*nX;
    nB = - diag(nx) + eye(n);
    nF=nAUX+nB;
    nF=(nF+nF')/2; % force symmetry
    [U,D]=eig(nF);
    lam=diag(D);
    nFinv=U*diag(1./lam)*U';
    ngrad=2*(diag(nFinv*nAUX)-nx);
    nres= norm(ngrad);

    [U,D]=eig(H);
    % sprintf("nbound:%f, nres:%f, Hessmineig:%f, min(nGamma):%f", nbound, nres,min(diag(D)),min(nGamma))

    if nbound-bound>c1*alfa*dir'*grad
        judge=0;
        b=1;
        a=0;
    elseif abs(dir'*ngrad)>c2*abs(dir'*grad)
        judge=0;
        b=2;
        a=0;
    else
        judge=1;
    end
    %line search
   
    while judge==0
        alfa=(a+b)/2;
        nGamma=Gamma.*exp(alfa*dir);
        [nbound,nx,ininfo] = SDPT3_BQP_light(nX,C,s,A_data,b_data,nGamma);
        nX=ininfo.X;
        scaleC=diag(nGamma)*C*diag(nGamma);
        nAUX = scaleC.*nX;
        nB = - diag(nx) + eye(n);
        nF=nAUX+nB;
        nF=(nF+nF')/2; % force symmetry
        [U,D]=eig(nF);
        lam=diag(D);
        nFinv=U*diag(1./lam)*U';
        ngrad=2*(diag(nFinv*nAUX)-nx);
        nres= norm(ngrad);
        if nbound-bound>c1*alfa*dir'*grad
            b=alfa;
        elseif abs(dir'*ngrad)>c2*abs(dir'*grad)
            a=alfa;
        else
            break
        end
        if abs(b-a)<1e-3
            break
        end   
    end
    deltay=log(nGamma./Gamma);
    Gamma=nGamma;
    deltag=ngrad-grad;
    grad=ngrad;
    res=nres;
    Finv=nFinv;
    B=nB;
    AUX=nAUX;
    X=nX;
    bound=nbound;
    gap=bound-heurval;
    allGamma=[allGamma,Gamma];
    allres=[allres,res];
    allbound=[allbound,bound];
    k=k+1; 
end

if gap <= TOL
    info.exitflag=0;
elseif abs(res) <= TOL
    info.exitflag=1;
elseif difgap <= TOL
    info.exitflag=2;
elseif k>Numiterations
    info.exitflag=3;
elseif toc(t1)> timelimit
    info.exitflag=4;
else
    info.exitflag=5;
end

info.maxiteration = Numiterations;
info.tol = TOL;
info.iterations=k-1;
info.gap=gap;
info.absres=abs(res);
info.difgap=difgap;
[optbound,optiteration]=min(allbound);
optGamma=allGamma(:,optiteration);
[bound1,~,~] = SDPT3_BQP_light(X0,C,s,A_data,b_data,ones(n,1));
if optbound>bound1
    optbound=bound1;
    optGamma=ones(n,1);
end
info.optbound=optbound;
info.optGamma=optGamma;
time=toc(t1);
info.time=time;