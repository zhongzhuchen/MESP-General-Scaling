%% obtain class properties
C=obj.C;
Cinv=obj.C_comp;
ldetC=obj.ldetC;
n=obj.size;
A_data=obj.A;
b_data=obj.b;
A_comp_data = -A_data;
b_comp_data = b_data-A_data*ones(n,1);
[m,~] = size(obj.A);
t1=tic;

if n>200
    TOL= 10^(-4);
    Numiterations=20; 
else
    TOL= 10^(-6);
    Numiterations=100; 
end

%% calculate the better lower bound among C and Cinv if C is invertible
heurval = obj.obtain_lb(s);

Y0=diag((n-s)/n*ones(n,1));
% Initialize Gamma
Gamma=GammaInit;

difgap=1;
k=1;

c1=1e-4;
c2=0.9;

%% calculate the gradient of linx bound with respect to Gamma
[bound,x,ininfo] = SDPT3_BQP_comp_light(Y0,Cinv,s,A_data,b_data,ldetC,Gamma);

Y=ininfo.Y;
y=ones(n,1)-x;
scaleCinv=diag(Gamma)*Cinv*diag(Gamma);
AUX = scaleCinv.*Y;
B = eye(n)-diag(y);
%Compute F(gamma,x)
F= AUX+B;
F=(F+F')/2; % force symmetry
%Compute inv(F(gamma,x))
[U,D]=eig(F);
lam=diag(D);
Finv=U*diag(1./lam)*U';
%Compute the residual res
gap=bound-heurval;

grad=2*(diag(Finv*AUX)+-y);
res= norm(grad);

allGamma=Gamma;
allres=res;
allbound=bound;

H=eye(n); % initialize the inverse Hessian approximation

%% we use the optimal solution of every last linx bound as the initial point 
% for solving the next linx bound, trick for accelarating optimization
nY=Y;

%% loop
while(k<=Numiterations && gap > TOL && abs(res) > TOL && difgap > TOL)
    sprintf('iteration: %d, res: %f',k,res)
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
    [nbound,nx,ininfo] = SDPT3_BQP_comp_light(nY,Cinv,s,A_data,b_data,ldetC,nGamma);
    nY=ininfo.Y;
    scaleCinv=diag(nGamma)*Cinv*diag(nGamma);
    nAUX = scaleCinv.*nY;
    nB = diag(nx);
    nF=nAUX+nB;
    nF=(nF+nF')/2; % force symmetry
    [U,D]=eig(nF);
    lam=diag(D);
    nFinv=U*diag(1./lam)*U';
    ngrad=2*(diag(nFinv*nAUX)+nx-ones(n,1));
    nres= norm(ngrad);

    [U,D]=eig(H);
    % sprintf("nbound:%f, nres:%f, Hessmineig:%f, min(nGamma):%f", nbound, nres,min(diag(D)),min(nGamma))

    if nbound-bound>c1*alfa*dir'*grad
        judge=0;
    elseif abs(dir'*ngrad)>c2*abs(dir'*grad)
        judge=0;
    else
        judge=1;
    end
    %line search
    b=1;
    a=0;
    while judge==0
        alfa=(a+b)/2;
        nGamma=Gamma.*exp(alfa*dir);
        [nbound,nx,ininfo] = SDPT3_BQP_comp_light(nY,Cinv,s,A_data,b_data,ldetC,nGamma);
        nY=ininfo.Y;
        scaleCinv=diag(nGamma)*Cinv*diag(nGamma);
        nAUX = scaleCinv.*nY;
        nB = diag(nx);
        nF=nAUX+nB;
        nF=(nF+nF')/2; % force symmetry
        [U,D]=eig(nF);
        lam=diag(D);
        nFinv=U*diag(1./lam)*U';
        ngrad=2*(diag(nFinv*nAUX)+nx-ones(n,1));
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
    Y=nY;
    bound=nbound;
    gap=bound-heurval;
    allGamma=[allGamma,Gamma];
    allres=[allres,res];
    allbound=[allbound,bound];
    k=k+1; 
end
info.iterations=k-1;
info.gap=gap;
info.absres=abs(res);
info.difgap=difgap;
[optbound,optiteration]=min(allbound);
optGamma=allGamma(:,optiteration);
[bound1,~,~] = SDPT3_BQP_comp_light(Y0,Cinv,s,A_data,b_data,ldetC,ones(n,1));
if optbound>bound1
    optbound=bound1;
    optGamma=ones(n,1);
end
info.optbound=optbound;
info.optGamma=optGamma;
time=toc(t1);
info.time=time;