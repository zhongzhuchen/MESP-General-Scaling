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

%% obtain the lower bound
heurval = obj.obtain_lb(s);

X0=diag(s/n*ones(n,1));

% Initialize gamma
csortoriginal=sort(diag(C),'descend');
gamma=1/csortoriginal(s); % initial scale factor if none provided
power=1.5;
gamma=gamma^power;

difgap=1;
k=1;

c1=1e-4;
c2=0.9;
% timelimit = 350;

%solve the linx ralaxation for gamma and obtain x
[bound,x,ininfo] = SDPT3_BQP_light(X0,C,s,A_data,b_data,sqrt(gamma)*ones(n,1));
X=ininfo.X;
AUX = C.*X;
%Compute F(gamma,x)
F=gamma*AUX - diag(x) + eye(n);
F=(F+F')/2; % force symmetry
%Compute inv(F(gamma,x))
[U,D]=eig(F);
lam=diag(D);
Finv=U*diag(1./lam)*U';
%Compute the residual res
gap=bound-heurval;
res= n - s - trace(Finv*diag(ones(n,1)-x)); 
allgamma=gamma;
allres=res;
allbound=bound;

nX=X;
while(k<=Numiterations && gap > TOL && abs(res) > TOL && difgap > TOL && toc(t1)<= timelimit) 
    sprintf('iteration: %d, res: %f',k,res);
    if k>1
        difgap=abs(allbound(k)-allbound(k-1));
    end
    %Compute MAT= H_G(gamma:=exp(psi))
    Mat = gamma*(ones(n,1)-x)'*diag(Finv*AUX*Finv);
    %Compute the search direction for Newton's method (dir=-res/Mat)
    dir=-res/Mat;
    edir=exp(dir);
    if abs(edir)==Inf
        edir=sign(edir)*1e10;
    end
    dir=log(edir);
    %check if alfa=1 satisfies the Strong Wolfe Conditions
    alfa=1;
    ngamma=gamma*exp(alfa*dir);
    [nbound,nx,ininfo] = SDPT3_BQP_light(nX,C,s,A_data,b_data,sqrt(ngamma)*ones(n,1));
    nX=ininfo.X;
    nAUX = C.*nX;
    nF=ngamma*nAUX - diag(nx) + eye(n);
    nF=(nF+nF')/2; % force symmetry
    [U,D]=eig(nF);
    lam=diag(D);
    nFinv=U*diag(1./lam)*U';
    nres= n - s - trace(nFinv*diag(ones(n,1)-nx)); 
    if nbound-bound>c1*alfa*dir*res
        judge=0;
    elseif abs(dir*nres)>c2*abs(dir*res)
        judge=0;
    else
        judge=1;
    end
    %line search
    b=1;
    a=0;
    while judge==0
        alfa=(a+b)/2;
        ngamma=gamma*exp(alfa*dir);
        [nbound,nx,ininfo] = SDPT3_BQP_light(nX,C,s,A_data,b_data,sqrt(ngamma)*ones(n,1));
        nX=ininfo.X;
        nAUX = C.*nX;
        nF=ngamma*nAUX - diag(nx) + eye(n);
        nF=(nF+nF')/2; % force symmetry
        [U,D]=eig(nF);
        lam=diag(D);
        nFinv=U*diag(1./lam)*U';
        nres= n - s - trace(nFinv*diag(ones(n,1)-nx)); 
        if nbound-bound>c1*alfa*dir*res
            b=alfa;
        elseif abs(dir*nres)>c2*abs(dir*res)
            a=alfa;
        else
            break
        end
        if abs(b-a)<1e-3
            break
        end   
    end
    gamma=ngamma;
    res=nres;
    Finv=nFinv;
    AUX=nAUX;
    x=nx;
    bound=nbound;
    gap=bound-heurval;
    allgamma=[allgamma,gamma];
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
optgamma=allgamma(optiteration);
[bound1,~,~] = SDPT3_BQP_light(X0,C,s,A_data,b_data,ones(n,1));
if optbound>bound1
    optbound=bound1;
    optgamma=1;
end
info.optbound=optbound;
info.optgamma=optgamma;
time=toc(t1);
info.time=time;
