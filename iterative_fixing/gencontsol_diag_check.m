function [xind] = gencontsol_diag_check(C,s,A,b,fix)
n=size(C,1);
f=log(diag(C));
vfix0 = fix.fixto0list;
vfix1 = fix.fixto1list;
IM = eye(n);
Aeq=[ones(1,n); IM(vfix1, :); IM(vfix0, :)];
beq=[s;ones(length(vfix1),1);zeros(length(vfix0),1)];
lb=zeros(n,1);
ub=ones(n,1);
OPTIONS.Display='none';
[x,fval,exitflag,OUTPUT] = intlinprog(f,1:n,A,b,Aeq,beq,lb,ub,[],OPTIONS);
if exitflag>=1
    xind=transpose(1:n);
    xind=xind(x>0.5);
else
    error('The MESP instance is infeasible or the solver does not find any feasible solution.')
end

% xType=2*ones(n,1);
% lb=zeros(n,1);
% ub=ones(n,1);
% x0=[];
% options = knitro_options('algorithm',3,...  % active-set/simplex algorithm
%                          'outlev',0);       % iteration display
% 
% [x,heurval,exitflag,~] = ...
%     knitro_milp(f,xType,A,b,Aeq,beq,lb,ub,x0,[],options);
% 
% if exitflag>-199 || -400 >= exitflag >= -406
%     xind=transpose(1:n);
%     xind=xind(x>0.5);
% elsei
%     error('The MESP instance is infeasible or the solver does not find any feasible solution.')
% end
end

