function [fval,dx,info] = DDFact_obj_Knitro_prescale(x,s,F,Fsquare,logGamma)
% light version DDFact_obj method of MESP by eliminating all unnessary
% calculations for saving time
% scale F with Gamma
n=length(x);
% column numbers of F
d=length(Fsquare(:,:,1));

X=zeros(d);
for i=1:n
    if x(i)==0
    else
        X=X+x(i)*Fsquare(:,:,i);
    end
end

% calculate k and corresponding mid_value as shown in Nikolov's paper
X=1/2*(X+X');
[U,D]=eig(X);
D=diag(D);
[sort_D, sort_idx]=sort(D, 'descend');
[k,mid_val]=find_k(sort_D,s);

if mid_val<=-n*1e-8 % for numerical stability
    sprintf('k=%d, mid_val=%f, sort_D(s)=%f, sum(x)=%f, rank(X)=%d', k, mid_val, sort_D(s), sum(x), rank(X))
    error('Something went wrong with calculating X or C might be a zero matrix.');
end

% construct the eigenvalue for the feasible solution to the dual problem, i.e., DFact
eigDual=zeros(d,1);

% for i=1:d
%     if(D(i)>mid_val)
%         eigDual(i)=1/D(i);
%     else
%         if mid_val>0
%             eigDual(i)=1/mid_val;
%         else
%             eigDual(i)=1/sort_D(k);
%         end
%     end
% end

% vectorized code of the above commented code
ind1=sort_idx(1:k);
non_ind1= setdiff(1:d, ind1);
eigDual(ind1)=1./D(ind1);
if mid_val>0
    eigDual(non_ind1)=1/mid_val;
else
    eigDual(non_ind1)=1/sort_D(k);
end

%% calculate the gradient dx
% dX=U*diag(eigDual)*U';
% dx=diag(F*dX*F');

% dX=U*diag(eigDual)*U';
% dx=diag(F*dX*F');

% vectorized code of the above commented code
K1=F*U;
K2=K1.*(eigDual');
dx1=sum(K2.*K1,2);
dx=dx1-logGamma;
sort_eigDual=sort(eigDual);
fval=-sum(log(sort_eigDual(1:s)))-sum(x.*logGamma);

fval=-fval;
dx=-dx;

info.cache = -s;
end
