%% obtain problem data
C=obj.C;
n = obj.size;
A_data=obj.A;
b_data=obj.b;
[m,~] = size(obj.A);
F=obj.F;
Fsquare=obj.Fsquare;
info=struct;
t1=tic;
%% setting different numbers of iterations for different problem size
if n>200
    TOL= 10^(-6);
    Numiterations=20; 
else
    TOL= 10^(-10);
    Numiterations=200; 
end

%% obtain the MESP optimal value lower bound
heurval = obj.obtain_lb(s);
x0=s/n*ones(n,1);
%% Initializes Gamma and hyperparameter
Gamma=GammaInit;

difgap=1;
k=1;
c1=1e-4;
c2=0.9;
% timelimit = 350;

%% calculate the gradient of fact bound with respect to Gamma
[bound,x,~]=Knitro_DDFact_light(x0,C,s,F,Fsquare,A_data,b_data,Gamma);

Fx=diag(sqrt(x))*F;
Fsquarex=Fsquare.*reshape(x,1,1,n);

[~,dGamma,~] = DDFact_obj_auxiliary(Gamma,s,Fx,Fsquarex);

gap=bound-heurval;
grad=Gamma.*dGamma-x;
res= norm(grad);

allGamma=Gamma;
allres=res;
allbound=bound;

H=eye(n); % initialize the inverse Hessian approximation
% sprintf('k: %d, gap: %f, abs(res): %f, difgap: %f',k, gap, abs(res), difgap)

%% we use the optimal solution of every last linx bound as the initial point 
% for solving the next linx bound, trick for accelarating optimization
nx=x;

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
    [nbound,nx,~]=Knitro_DDFact_light(nx,C,s,F,Fsquare,A_data,b_data,nGamma);
    nFx=diag(sqrt(nx))*F;
    nFsquarex=Fsquare.*reshape(nx,1,1,n);
    [~,dnGamma,~] = DDFact_obj_auxiliary(nGamma,s,nFx,nFsquarex);

    ngrad=nGamma.*dnGamma-nx;
    nres= norm(ngrad);
    % [U,D]=eig(H);
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
        [nbound,nx,~]=Knitro_DDFact_light(nx,C,s,F,Fsquare,A_data,b_data,nGamma);
        nFx=diag(sqrt(nx))*F;
        nFsquarex=Fsquare.*reshape(nx,1,1,n);
        [~,dnGamma,~] = DDFact_obj_auxiliary(nGamma,s,nFx,nFsquarex);
        ngrad=nGamma.*dnGamma-nx;
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
    x=nx;
    bound=nbound;
    gap=bound-heurval;
    allGamma=[allGamma,Gamma];
    allres=[allres,res];
    allbound=[allbound,bound];
    k=k+1; 

    %sprintf('k: %d, gap: %f, abs(res): %f, difgap: %f',k, gap, abs(res), difgap)
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
[bound1,~,~]=Knitro_DDFact_light(x0,C,s,F,Fsquare,A_data,b_data,ones(n,1));
if optbound>bound1
    optbound=bound1;
    optGamma=ones(n,1);
end
info.optbound=optbound;
info.optGamma=optGamma;
time=toc(t1);
info.time=time;