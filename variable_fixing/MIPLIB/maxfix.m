function [fix0vec, fix1vec, info] = maxfix(problem, optval)
%{
This is a function exploring the maximum power of variable fixing for some
upper bound of MESP. In other words, given an index i, this function try to
detect if there is any dual feasible solution that can fix x_i either to 0
or 1 under the fixing pattern of the paper.

info: a struct containing output information

%}
% collect problem data
f = problem.f;
A = problem.Aineq;
b = problem.bineq;
Aeq = problem.Aeq;
beq = problem.beq;
lb = problem.lb;
ub = problem.ub;
intcon = problem.intcon;
% detect binary variables
binarycon = find((lb==0) & (ub ==1));
binarycon = intersect(binarycon,intcon);

n = length(f);
[m1,~] = size(A);
[m2,~] = size(Aeq);
% array for storing fixed variables
fix0vec = [];
fix1vec = [];
info = struct;

% record total time used
TStart=tic;
tStart=cputime;
% options for Gurobi
options.OutputFlag = 0;

fix0tol = [];
fix1tol = [];
binarycon = binarycon';
for ind =  binarycon
    sprintf("Check index: %d, elapse time: %f.", ind, toc(TStart))
    lb=zeros(n,1);
    ub=ones(n,1);
    lb(ind) = 1-1e-12;
    [x,fval,exitflag,output,lambda, result] = linprog(f,A,b,Aeq,beq,lb,ub,[],options);
    % check feasibility / fixtol
    if exitflag == -2 || (exitflag == 1 && fval >= optval+1e-10)
        fix0vec(end+1) = ind; 
    end
    if exitflag == -2
        fix0tol(end+1) = 0;
    elseif exitflag == -3
        fix0tol(end+1) = Inf;
    else
        fix0tol(end+1) = fval-optval;
    end
    
    % check if x_ind can be fixed to 1
    lb=zeros(n,1);
    ub=ones(n,1);
    ub(ind) = 1e-12;
    [x,fval,exitflag,output,lambda, result] = linprog(f,A,b,Aeq,beq,lb,ub,[],options);
    % check feasibility / fixtol
    if exitflag == -2 || (exitflag == 1 && fval >= optval+1e-10)
        fix1vec(end+1) = ind; 
    end
    if exitflag == -2
        fix1tol(end+1) = 0;
    elseif exitflag == -3
        fix1tol(end+1) = Inf;
    else
        fix1tol(end+1) = fval-optval;
    end
end
time=toc(TStart);
tEnd=cputime-tStart;

info.optval=optval;
info.fix0vec = fix0vec;
info.fix1vec = fix1vec;
info.fixnum = length(fix0vec)+length(fix1vec);
info.time = time;
info.cputime = tEnd;
info.fix0tol = fix0tol;
info.fix1tol = fix1tol;
end