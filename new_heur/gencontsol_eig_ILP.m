function [xind] = gencontsol_eig_ILP(C,s,A,b)
n=size(C,1);
try
    [v,e]=eig(C);
catch
    [v,e]=eigs(C,n);
end
v=transpose(v);
f=sum(v(:, (n-s+1):n).^2, 2);

Aeq=ones(1,n);
beq=s;
xType=2*ones(n,1);
lb=zeros(n,1);
ub=ones(n,1);
x0=[];
options = knitro_options('algorithm',3,...  % active-set/simplex algorithm
                         'outlev',0);       % iteration display

[x,heurval,exitflag,~] = ...
    knitro_milp(f,xType,A,b,Aeq,beq,lb,ub,x0,[],options);

if exitflag>-199 || -400 >= exitflag >= -406
    xind=transpose(1:n);
    xind=xind(x>0.5);
else
    error('The MESP instance is infeasible or the solver does not find any feasible solution.')
end
end

