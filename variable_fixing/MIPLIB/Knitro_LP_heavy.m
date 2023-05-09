function [info] = Knitro_LP_heavy(problem, optval, info_input)
%{
we try to fix with dual feasible solutions at each iteration where problem 
is a structure read from the MIPLIB mps file that can be solved by MATLAB 
'intlinprog', optval is the optimal value of the corresponding problem.

Example of problem structure:

          f: [28926×1 double]
      Aineq: [21915×28926 double]
      bineq: [21915×1 double]
        Aeq: [1379×28926 double]
        beq: [1379×1 double]
         lb: [28926×1 double]
         ub: [28926×1 double]
     intcon: [28926×1 double]
     solver: 'intlinprog'
    options: [1×1 optim.options.Intlinprog]

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
if m1 == 0
    A = [];
    b = [];
end
if m2 == 0
    Aeq = [];
    beq = [];
end

% array for storing fixed variables
fix0vec = [];
fix1vec = [];

itercount = 0;
% structure for storing results
info = struct;
% define nested objective function for Knitro which fix with dual feasible
% solution at every iteration
    function terminate = outfun(x,optimValues,state)
        if ~isinf(optimValues.fval)
            itercount = itercount +1;
            % sprintf('Iteration: %d', itercount)
            % obtain dual solution
            upsilon = optimValues.lambda.lower;
            nu = optimValues.lambda.upper;
            pi = optimValues.lambda.ineqlin;
            tau = optimValues.lambda.eqlin;
%             % check primal feasibility
%             if (norm(max(A*x-b,0))<1e-10*m1) &...
%                     (norm(Aeq*x-beq)<1e-10*m2) &...
%                     (norm(min(x-lb,0))<1e-10*n) &...
%                     (norm(max(x-ub,0))<1e-10*n)
%                 sprintf('At iteration: %d, we have primal feasible solution.', itercount)
%             end
            % check dual feasibility
            if m1 == 0 && m2 > 0
                dualresidual_eq = norm(f+Aeq'*tau+upsilon+nu);
                lowerbound = -tau'*beq-sum(nu);
            elseif m1 > 0 && m2 == 0
                dualresidual_eq = norm(f+A'*pi+upsilon+nu);
                lowerbound = -pi'*b-sum(nu);
            elseif m1 == 0 && m2 == 0
                dualresidual_eq = norm(f+upsilon+nu);
                lowerbound = -sum(nu);
            else
                 dualresidual_eq = norm(f+A'*pi+Aeq'*tau+upsilon+nu);
                 lowerbound = -pi'*b-tau'*beq-sum(nu);
            end
            if (norm(max(upsilon,0)) < 1e-10*n) &...
                    (norm(min(nu,0)) < 1e-10*n) &...
                    (norm(min(pi,0)) < 1e-10*m1) &...
                    (dualresidual_eq<1e-10*n)
                sprintf('At iteration: %d, we have dual feasible solution.', itercount)
                lowerbound = min(lowerbound, optval);
                % check variables that can be fixed to 0
                indicator = lowerbound - upsilon -optval - 1e-10 > 0;
                ind = find(indicator);
                ind = intersect(binarycon, ind);
                old0len = length(fix0vec);
                fix0vec = union(fix0vec,ind);
                % check variables that can be fixed to 1
                indicator = lowerbound + nu -optval - 1e-10 > 0;
                ind = find(indicator);
                ind = intersect(binarycon, ind);
                old1len = length(fix1vec);
                fix1vec = union(fix1vec,ind);
                sprintf('At iteration: %d, %d more variables are fixed.', itercount,...
                    length(fix0vec)+length(fix1vec)-old0len-old1len)
            end
        end
        terminate = false;
    end
% solve the linear programming problem with dualsimplex method
options = knitro_options('algorithm',1,... % active-set/simplex algorithm
               'outlev',0,...   % iteration display
               'presolve',0,... % deactive Knitro pre-solver
               'feastol', 1e-10,... % feasibility tolerance
               'opttol', 1e-10,...
               'bar_feasible', 0); 
extendedFeatures.OutputFcn = @ outfun;
% initial lambda value if provided
if exist('info_input', 'var')
    extendedFeatures.lambdaInitial.eqlin = info_input.tau;
    extendedFeatures.lambdaInitial.ineqlin = info_input.pi;
    extendedFeatures.lambdaInitial.lower = -info_input.upsilon;
    extendedFeatures.lambdaInitial.upper = info_input.nu;
end
[x,fval,exitflag,output,lambda] = knitro_lp (f, A, b, Aeq, beq, lb, ub,[], extendedFeatures, options);
% fix with the final solution
fix0vec_final = [];
fix1vec_final = [];
upsilon = lambda.lower;
nu = lambda.upper;
pi = lambda.ineqlin;
tau = lambda.eqlin;
if (0 >= exitflag && exitflag > -150) || (-300 >= exitflag &&exitflag >= -406)
    info.feasibility = 1;
    % check dual feasibility 
    if m1 == 0 && m2 > 0
        dualresidual_eq = norm(f+Aeq'*tau+upsilon+nu);
        lowerbound = -tau'*beq-sum(nu);
    elseif m1 > 0 && m2 == 0
        dualresidual_eq = norm(f+A'*pi+upsilon+nu);
        lowerbound = -pi'*b-sum(nu);
    elseif m1 == 0 && m2 == 0
        dualresidual_eq = norm(f+upsilon+nu);
        lowerbound = -sum(nu);
    else
         dualresidual_eq = norm(f+A'*pi+Aeq'*tau+upsilon+nu);
         lowerbound = -pi'*b-tau'*beq-sum(nu);
    end
    if (norm(max(upsilon,0)) < 1e-10*n) &...
            (norm(min(nu,0)) < 1e-10*n) &...
            (norm(min(pi,0)) < 1e-10*m1) &...
            (dualresidual_eq<1e-10*n)
        lowerbound = min(lowerbound, optval);
        % check variables that can be fixed to 0
        indicator = lowerbound - upsilon -optval - 1e-10 > 0;
        ind = find(indicator);
        ind = intersect(binarycon, ind);
        fix0vec_final = union(fix0vec_final,ind);
        % check variables that can be fixed to 1
        indicator = lowerbound + nu -optval - 1e-10 > 0;
        ind = find(indicator);
        ind = intersect(binarycon, ind);
        fix1vec_final = union(fix1vec_final,ind);
    end
end
info.exitflag = exitflag;
info.iterations = output.iterations;
info.fval = fval;
info.x = x;
info.integrality_gap = optval - fval;
info.fix0vec_prefinal = setdiff(fix0vec, fix0vec_final);
info.fix1vec_prefinal = setdiff(fix1vec, fix1vec_final);
info.fix0vec_final = fix0vec_final;
info.fix1vec_final = fix1vec_final;
info.fix0vec = fix0vec;
info.fix1vec = fix1vec;
info.fixvec = union(fix0vec, fix1vec);
info.upsilon = -upsilon;
info.nu = nu;
info.pi = pi;
info.tau= tau;
end

