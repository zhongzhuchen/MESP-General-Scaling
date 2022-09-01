%% obtain class properties
C=obj.C;
n = obj.size;
info=struct;

t1=tic;
if n>200
    TOL= 10^(-6);
    Numiterations=20; 
else
    TOL= 10^(-10);
    Numiterations=200; 
end

%Tolerance for residual and dualtity gap
%Maximum number of Newton steps to optimize psi=log(gamma)
[U,D]=eig(C);
lam=diag(D);

% calculate the better lower bound among C and Cinv if C is invertible
if min(lam)> 0
    shift=log(prod(lam));  % logdet C
    Cinv=U*diag(1./lam)*U';
    [~,heurval]=heur(C,n,s);        % HEURSITIC ON ORIGINAL
    [~,cheurval]=heur(Cinv,n,n-s); % HEURISTIC ON COMPLEMENT
    if cheurval+shift > heurval        % PICK THE BEST
        heurval=cheurval+shift;
    end
else
    [~,heurval]=heur(C,n,s);        % HEURSITIC ON ORIGINAL
end

x0=s/n*ones(n,1);
% Initialize Gamma
Gamma=GammaInit;

difgap=1;
k=1;

c1=1e-4;
c2=0.9;
%solve the linx ralaxation for gamma and obtain x
[bound,x,~]= obj.Knitro_Linx(x0,s,Gamma);
scaleC=diag(Gamma)*C;
AUX = scaleC*diag(x)*scaleC';
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

grad=diag(Finv*AUX)-x;
res= norm(grad);

allGamma=Gamma;
allres=res;
allbound=bound;

H=eye(n); % initialize the inverse Hessian approximation
while(k<=Numiterations && gap > TOL && abs(res) > TOL && difgap > TOL)
    sprintf('iteration: %d',k);
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
    if norm(edir)==Inf || isnan(norm(edir))
        edir=edir/norm(edir)*1e10;
    end
    dir=log(edir);
    %check if alfa=1 satisfies the Strong Wolfe Conditions
    alfa=1;
    nGamma=Gamma.*exp(alfa*dir);
    [nbound,nx,~]= obj.Knitro_Linx(x0,s,nGamma);
    scaleC=diag(nGamma)*C;
    nAUX = scaleC*diag(nx)*scaleC';
    nB = - diag(nx) + eye(n);
    nF=nAUX+nB;
    nF=(nF+nF')/2; % force symmetry
    [U,D]=eig(nF);
    lam=diag(D);
    nFinv=U*diag(1./lam)*U';
    ngrad=diag(nFinv*nAUX)-nx;
    nres= norm(ngrad);

    [U,D]=eig(H);
    sprintf("nbound:%f, nres:%f, Hessmineig:%f, min(nGamma):%f", nbound, nres,min(diag(D)),min(nGamma))

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
        [nbound,nx,~]= obj.Knitro_Linx(x0,s,nGamma);
        scaleC=diag(nGamma)*C;
        nAUX = scaleC*diag(nx)*scaleC';
        nB = - diag(nx) + eye(n);
        nF=nAUX+nB;
        nF=(nF+nF')/2; % force symmetry
        [U,D]=eig(nF);
        lam=diag(D);
        nFinv=U*diag(1./lam)*U';
        ngrad=diag(nFinv*nAUX)-nx;
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
    x=nx;
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
[bound1,~,~]= obj.Knitro_Linx(x0,s,ones(n,1));
if optbound>bound1
    optbound=bound1;
    optGamma=ones(n,1);
end
info.optbound=optbound;
info.optGamma=optGamma;
time=toc(t1);
info.time=time;