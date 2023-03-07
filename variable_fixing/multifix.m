function [fix0vec, fix1vec, info] = multifix(ubname, C, s, A, b, I, J, x0, Gamma, LB)
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
    I: set of indices aimed to fix to 0
    J: set of indices aimed to fix to 1
    x0: inital point
    Gamma: general scaling parameter
    LB: lower bound for the MESP instance

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
info = struct;
fix0vec = [];
fix1vec = [];
% set up constraints
lb=zeros(n,1);
ub=ones(n,1);
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

lb = zeros(n,1);
ub = ones(n,1);
lb(I) = 1/(length(I)+length(J));
ub(J) = 1-1/(length(I)+length(J));

if  strcmp(ubname, 'linx')
    [x,knitro_fval,exitflag,output,lambda,~] = knitro_nlp(obj_fn,x0,A,b,Aeq,beq,lb,ub,[],extendedFeatures,options);
elseif (strcmp(ubname, 'Fact')) || (strcmp(ubname, 'cFact'))
    [x,knitro_fval,exitflag,output,lambda,~] = knitro_nlp(obj_fn,x0,A,b,Aeq,beq,lb,ub,[],[],options);
end
time=toc(TStart);
tEnd=cputime-tStart;

if (0 >= exitflag && exitflag > -150) || (-300 >= exitflag &&exitflag >= -406)
    info.feasibility = 1;
    %% construct dual feasible solutions
    if strcmp(ubname, 'linx')
        [fval, dx, ininfo] = Linx_obj_Knitro(x,C,Gamma);
        dx = -dx; % Note the subroutine produce negative sign, converts back
        fval = -fval;
        dx2 = ininfo.dx2;
    elseif strcmp(ubname, 'Fact')
        [fval,dx,~] = DDFact_obj_Knitro(x,s,F,Fsquare,Gamma);
        dx = -dx;
        fval = -fval;
    elseif strcmp(ubname, 'cFact')
        [fval,dx,~] = DDFact_comp_obj_Knitro(x,s,F_comp,Fsquare_comp,ldetC,Gamma);
        dx=-dx;
        fval = -fval;
    end
    
    c1 = zeros(n,1);
    c1(I) = -1/(length(I)+length(J));
    c2 = ones(n,1);
    c2(J) = 1-1/(length(I)+length(J));
    f=[c1;c2;b;s];
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
        dualbound = dualgap+1/(length(I)+length(J))*(sum(info.dual_upsilon(I))+sum(info.dual_nu(J)))...
            +0.5*sum(dx2)-n/2+sum(x.*log(Gamma))+fval;
    elseif strcmp(ubname, 'Fact')
        dualbound = dualgap+1/(length(I)+length(J))*(sum(info.dual_upsilon(I))+sum(info.dual_nu(J)))...
            -s+sum(x.*log(Gamma))+fval;
    elseif strcmp(ubname, 'cFact')
        dualbound = dualgap+1/(length(I)+length(J))*(sum(info.dual_upsilon(I))+sum(info.dual_nu(J)))...
            -(n-s)-sum(dx)+sum((ones(n,1)-x).*log(Gamma))+fval;
    end
else
    info.feasibility = 0;
    info.dual_upsilon = zeros(n,1);
    info.dual_nu = zeros(n,1);
    dualbound = Inf;
end

info.integrality_gap= dualbound-LB;
%% fix variables
%% fixing variables by the fixing logic in the factorization paper by Chen etal.
info.fixnum=0;
info.fixnum_to0=0;
info.fixto0list=[];
info.fixnum_to1=0;
info.fixto1list=[];

% if the integrality gap between the upper bound and the lower bound is
% large, the fixing power might be hidden due the weak lower bound

if info.integrality_gap>1e-6
    info.solved=0;
else
    info.solved=1;
end
for i=1:n
    if info.integrality_gap<info.dual_upsilon(i)-1e-10 % fix to zero
        info.fixnum=info.fixnum+1;
        info.fixnum_to0=info.fixnum_to0+1;
        info.fixto0list(end+1)=i;
    elseif info.integrality_gap<info.dual_nu(i)-1e-10 % fix to one
        info.fixnum=info.fixnum+1;
        info.fixnum_to1=info.fixnum_to1+1;
        info.fixto1list(end+1)=i;
    end
end

% % info.fixto0array = zeros(n,1);
% % info.fixto0array(info.fixto0list) = 1;
% % info.fixto1array = zeros(n,1);
% % info.fixto1array(info.fixto1list) = 1;
% 
% % implement conflict matrices here (not the most efficient way):
% %{
% We store four n-by-n matrices here, where if the (i,j) element of:
% 1. the 1st matrix is one, then x_i and x_j cannot be one/ one simultaneously
% 2. the 2nd matrix is one, then x_i and x_j cannot be one/ zero simultaneously
% 3. the 3rd matrix is one, then x_i and x_j cannot be zero/ one simultaneously
% 4. the 4th matrix is one, then x_i and x_j cannot be zero/ zero simultaneously
% %}
% A11 = (info.integrality_gap < info.dual_upsilon + info.dual_upsilon' -1e-10);
% A10 = (info.integrality_gap < info.dual_upsilon + info.dual_nu' -1e-10);
% A01 = (info.integrality_gap < info.dual_nu + info.dual_upsilon' -1e-10);
% A00 = (info.integrality_gap < info.dual_nu + info.dual_nu' -1e-10);
% 
% % leverage the above conflict matrices to induce more variable fixing
% fixto0list = info.fixto0list;
% fixto1list = info.fixto1list;
% while true
%     oldlength0 = length(fixto0list);
%     oldlength1 = length(fixto1list);
%     % fix more variables to 0
%     A11sub_colsum = sum(A11(:, fixto1list), 2) > 0;
%     fixto0list = union(fixto0list, find(A11sub_colsum));
%     A10sub_colsum = sum(A10(:, fixto0list), 2) > 0;
%     fixto0list = union(fixto0list, find(A10sub_colsum));
%     % fix more variables to 0
%     A01sub_colsum = sum(A01(:, fixto1list), 2) > 0;
%     fixto1list = union(fixto1list, find(A01sub_colsum));
%     A00sub_colsum = sum(A00(:, fixto0list), 2) > 0;
%     fixto1list = union(fixto1list, find(A00sub_colsum));
%     if oldlength0 == length(fixto0list) && oldlength1 == length(fixto1list)
%         break;
%     end
% end
% 
% fixto0list = unique(fixto0list);
% fixto1list = unique(fixto1list);

% info.fix1num = length(fixto1list);
% info.fixto1list = fixto1list;
% fix1vec = fixto1list;
% info.fix0num = length(fixto0list);
% info.fixto0list = fixto0list;
% fix0vec = fixto0list;
fix1vec = info.fixto1list;
fix0vec = info.fixto0list;
info.fix0num = length(info.fixto0list);
info.fix1num = length(info.fixto1list);

info.fixnum = info.fix0num + info.fix1num;
info.time = time;
info.cputime = tEnd;

end