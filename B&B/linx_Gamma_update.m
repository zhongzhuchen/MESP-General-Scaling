c1=1e-4;
c2=0.9;
A_data = double.empty(0,n);
b_data = double.empty(0,1);
nx=xval;
grad=diag(Finv*AUX)-xval;
dir=-grad;
edir=exp(dir);
if norm(edir)==Inf || norm(edir)==0 || isnan(norm(edir))
    if norm(dir)==Inf || isnan(norm(dir))
        dir=-grad;
    else
        dir=dir/norm(dir)*10;
    end   
end
alfa=1;
nGamma=Gamma.*exp(alfa*dir);
[nbound,nx,~] = Knitro_Linx_light(nx,C,s,A_data,b_data,nGamma);
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
    [nbound,nx,~] = Knitro_Linx_light(nx,C,s,A_data,b_data,nGamma);
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
newGamma=nGamma; % update Gamma