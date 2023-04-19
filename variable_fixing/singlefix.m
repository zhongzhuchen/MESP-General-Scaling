function [fix0vec, fix1vec, info] = singlefix(ubname, testind, C, s, A, b, x0, Gamma, LB, xhat)
%{
Input:
    ub: the name of the upper bound, choose from
        - 'linx'
        - 'Fact'
        - 'cFact'
    testind: the set of indices for each if there is any dual feasible
        solution that can fix it
    C, s, A, b: problem data where C is the covariance matrix, s is the
        cardinality of the subset to choose, A and b are parameters for
        linear constraints
    x0: initial point
    Gamma: general scaling parameter
    LB: lower bound for the MESP instance
    xhat: the given fixed xhat, which will give a corresponding Theta, note
        that xhat does not have to be optimal for the primal problem
    
Output:
    fix0vec: the collection of indices that can be fixed to 0
    fix1vec: the collection of indices that can be fixed to 1
    info: a struct containing additional information
%}
n = length(C);
[m, ~] = size(A);
if ~exist('x0', 'var')
    x0 = s/n*ones(n,1);
end
if ~exist('Gamma', 'var')
    Gamma = ones(n,1);
end
if ~exist('LB', 'var')
    MESPEX1 = MESP(C, A, b);
    LB = MESPEX1.obtain_lb(s);
end
% if xhat is not given, construct Theta by an optimal solution to the
% primal problem
if ~exist('xhat', 'var')
    if strcmp(ubname, 'linx')
        [~,x,~] = Knitro_Linx_light(x0,C,s,A,b,Gamma);
        [fval, dx, ininfo] = Linx_obj_Knitro(x,C,Gamma);
        dx = -dx; % Note the subroutine produce negative sign, converts back
        fval = -fval;
        dx2 = ininfo.dx2;
    elseif strcmp(ubname, 'Fact')
        [F,Fsquare,ldetC] = gen_data(C,0);
        [~,x,info]=Knitro_DDFact_light(x0,C,s,F,Fsquare,A,b,Gamma);
        [fval,dx,~] = DDFact_obj_Knitro(x,s,F,Fsquare,Gamma);
        dx = -dx;
        fval = -fval;
    elseif strcmp(ubname, 'cFact')
        [F_comp,Fsquare_comp,ldetC] = gen_data(C,1);
        [~,x,info] = Knitro_DDFact_comp_light(x0,C,s,F_comp,Fsquare_comp,ldetC,A,b,Gamma);
        [fval,dx,~] = DDFact_comp_obj_Knitro(x,s,F_comp,Fsquare_comp,ldetC,Gamma);
        dx=-dx;
        fval = -fval;
    end
else
    if strcmp(ubname, 'linx')
        [fval, dx, ininfo] = Linx_obj_Knitro(xhat,C,Gamma);
        dx = -dx; % Note the subroutine produce negative sign, converts back
        fval = -fval;
        dx2 = ininfo.dx2;
    elseif strcmp(ubname, 'Fact')
        [F,Fsquare,ldetC] = gen_data(C,0);
        [fval,dx,~] = DDFact_obj_Knitro(xhat,s,F,Fsquare,Gamma);
        dx = -dx;
        fval = -fval;
    elseif strcmp(ubname, 'cFact')
        [F_comp,Fsquare_comp,ldetC] = gen_data(C,1);
        [fval,dx,~] = DDFact_comp_obj_Knitro(xhat,s,F_comp,Fsquare_comp,ldetC,Gamma);
        dx=-dx;
        fval = -fval;
    end
end

info = struct;
fix0vec = [];
fix1vec = [];

% the linear reduced problem for all upper bounds are the same
f=[zeros(n,1);ones(n,1);b;s];
Aeq=[-eye(n),eye(n),A',ones(n,1)];
beq=dx;
lb=[zeros(2*n+m,1);-inf];
ub=Inf(2*n+m+1,1);
options = knitro_options('algorithm',3,...  % active-set/simplex algorithm
                         'outlev',0);       % iteration display
[xlp, dualgap, exitflag, ~] = knitro_lp (f, [], [], Aeq, beq, lb, ub, [], [], options);

info.dual_upsilon = xlp(1:n);
info.dual_nu = xlp((n+1):2*n);
if strcmp(ubname, 'linx')
    dualbound = dualgap+0.5*sum(dx2)-n/2+sum(x.*log(Gamma))+fval;
elseif strcmp(ubname, 'Fact')
    dualbound = dualgap-s+sum(x.*log(Gamma))+fval;
elseif strcmp(ubname, 'cFact')
    dualbound = dualgap-(n-s)-sum(dx)+sum((ones(n,1)-x).*log(Gamma))+fval;
end
info.integrality_gap= dualbound-LB;

% record total time used
% Todo: improving the way for obtaining a better single dual feasible
% solution for variable fixing
TStart=tic;
tStart=cputime;

for ind=testind
    % check if x_ind can be fixed to 0
    if dualbound - xlp(ind) <= LB -1e-10
        fix0vec(end+1) = ind;
    end
    % check if x_ind can be fixed to 1
    if dualbound - xlp(n+ind) <= LB -1e-10
        fix1vec(end+1) = ind;
    end
end

info.fix0num = length(fix0vec);
info.fix1num = length(fix1vec);
info.fixto0list = fix0vec;
info.fixto1list = fix1vec;


time=toc(TStart);
tEnd=cputime-tStart;

info.fixnum = length(fix0vec)+length(fix1vec);
info.time = time;
info.cputime = tEnd;

end