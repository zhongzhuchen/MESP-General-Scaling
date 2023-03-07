function [fix0vec, fix1vec, info] = maxfix(ubname, testind, C, s, A, b, x0, Gamma, LB)
%{
This is a function exploring the maximum power of variable fixing for some
upper bound of MESP. In other words, given an index i, this function try to
detect if there is any dual feasible solution that can fix x_i either to 0
or 1 under the fixing pattern of the paper.

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
    x0: inital point
    Gamma: general scaling parameter
    LB: lower bound for the MESP instance

Output:
    fix0vec: the collection of indices that can be fixed to 0
    fix1vec: the collection of indices that can be fixed to 1
    info: a struct containing additional information
%}
n = length(C);
[m,~] = size(A);
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
info = struct;
fix0vec = [];
fix1vec = [];
% set up constraints
Aeq=ones(1,n);
beq=s;


% set up objective function based on ub
if strcmp(ubname, 'linx')
    obj_fn =  @(x) Linx_obj_Knitro(x,C,Gamma);
    extendedFeatures.HessFcn = @(x,lambda) Linx_hessfun(x,lambda,diag(Gamma)*C);
    options = knitro_options('algorithm', 0, 'convex', 1, 'derivcheck', 0, 'outlev', 0 , 'gradopt', 1, ...
                         'hessopt', 1, 'maxit', 1000, 'xtol', 1e-15, 'derivcheck_tol',1e-5,...
                         'feastol', 1e-10, 'opttol', 1e-10, 'bar_feasible',1,...
                         'bar_maxcrossit', 10);
elseif strcmp(ubname, 'Fact')
    [F,Fsquare,ldetC] = gen_data(C,0);
    obj_fn =  @(x) DDFact_obj_Knitro(x,s,F,Fsquare,Gamma);
    options = knitro_options('algorithm', 0, 'convex', 1, 'derivcheck', 0, 'outlev', 0 , 'gradopt', 1, ...
                         'hessopt', 2, 'maxit', 1000, 'xtol', 1e-15, ...
                         'feastol', 1e-10, 'opttol', 1e-10, 'bar_feasible',1,...
                         'bar_maxcrossit', 10);
elseif strcmp(ubname, 'cFact')
    [F_comp,Fsquare_comp,ldetC] = gen_data(C,1);
    obj_fn =  @(x) DDFact_comp_obj_Knitro(x,s,F_comp,Fsquare_comp,ldetC,Gamma);
    options = knitro_options('algorithm', 0, 'convex', 1, 'derivcheck', 0, 'outlev', 0 , 'gradopt', 1, ...
                         'hessopt', 2, 'maxit', 1000, 'xtol', 1e-15, ...
                         'feastol', 1e-10, 'opttol', 1e-10, 'bar_feasible',1,...
                         'bar_maxcrossit', 10);
end

% record total time used
TStart=tic;
tStart=cputime;

fix0tol = [];
fix1tol = [];
I = eye(n);
for ind=testind
    lb=zeros(n,1);
    ub=ones(n,1);
    lb(ind) = 1-1e-10;
    if  strcmp(ubname, 'linx')
        [x,knitro_fval,exitflag,output,lambda,~] = knitro_nlp(obj_fn,x0,A,b,Aeq,beq,lb,ub,[],extendedFeatures,options);
    elseif (strcmp(ubname, 'Fact'))
        [x,knitro_fval,exitflag,output,lambda,~] = knitro_nlp(obj_fn,x0,A,b,Aeq,beq,lb,ub,[],[],options);
    elseif(strcmp(ubname, 'cFact'))
        [x,knitro_fval,exitflag,output,lambda,~] = knitro_nlp(obj_fn,x0,A,b,Aeq,beq,lb,ub,[],[],options);
    end
    
    % check feasibility / fixtol
    if ~((0 >= exitflag && exitflag > -150) || (-300 >= exitflag &&exitflag >= -406)) || -knitro_fval <= LB
        fix0vec(end+1) = ind; 
    end
    fix0tol(end+1) = [-knitro_fval - LB];
    % check if x_ind can be fixed to 1
    lb=zeros(n,1);
    ub=ones(n,1);
    ub(ind) = 1e-10;
    if  strcmp(ubname, 'linx')
        [x,knitro_fval,exitflag,output,lambda,~] = knitro_nlp(obj_fn,x0,A,b,Aeq,beq,lb,ub,[],extendedFeatures,options);
    elseif (strcmp(ubname, 'Fact'))
        [x,knitro_fval,exitflag,output,lambda,~] = knitro_nlp(obj_fn,x0,A,b,Aeq,beq,lb,ub,[],[],options);
    elseif(strcmp(ubname, 'cFact'))
        [x,knitro_fval,exitflag,output,lambda,~] = knitro_nlp(obj_fn,x0,A,b,Aeq,beq,lb,ub,[],[],options);
    end
    
    % check feasibility / fixtol
    if ~((0 >= exitflag && exitflag > -150) || (-300 >= exitflag &&exitflag >= -406)) || -knitro_fval <= LB
        fix1vec(end+1) = ind; 
    end
    fix1tol(end+1) = [-knitro_fval - LB];
end
time=toc(TStart);
tEnd=cputime-tStart;

info.fixto0list = fix0vec;
info.fixto1list = fix1vec;
info.fix0num = length(fix0vec);
info.fix1num = length(fix1vec);
info.fixnum = length(fix0vec)+length(fix1vec);
info.time = time;
info.cputime = tEnd;

info.fix0tol = fix0tol;
info.fix1tol = fix1tol;
end