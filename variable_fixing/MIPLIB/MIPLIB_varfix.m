% problem = mpsread('50v-10.mps'); optval = 1;
% problem = mpsread('academictimetablesmall.mps'); optval = 0;
% problem = mpsread('air05.mps'); optval = 26374;
% problem = mpsread('bppc4-08.mps'); optval = 53;
% problem = mpsread('chromaticindex512-7.mps'); optval = 4;
% problem = mpsread('chromaticindex32-8.mps'); optval = 4;
% problem = mpsread('2club200v15p5scn.mps'); optval = -70;

% f = problem.f;
% A = problem.Aineq;
% b = problem.bineq;
% Aeq = problem.Aeq;
% beq = problem.beq;
% lb = problem.lb;
% ub = problem.ub;
% 
% % options.IterationLimit = 10000;
% options.Presolve = 0;
% options.Method = 2;
% options.OutputFlag = 1;
% [x,fval,exitflag,output,lambda, result] = linprog(f,A,b,Aeq,beq,lb,ub,[],options);
% 
% info = struct;
% 
% n = length(f);
% [m1,~] = size(A);
% [m2,~] = size(Aeq);
% 
% fix0vec_final = [];
% fix1vec_final = [];
% 
% upsilon = lambda.lower;
% nu = lambda.upper;
% pi = lambda.ineqlin;
% tau = lambda.eqlin;
% % check dual feasibility 
% if (norm(min(upsilon,0)) < 1e-10*n) &...
%         (norm(min(nu,0)) < 1e-10*n) &...
%         (norm(min(pi,0)) < 1e-10*m1) &...
%         (norm(f+A'*pi+Aeq'*tau-upsilon+nu)<1e-10*n)
%     lowerbound = -pi'*b-tau'*beq-sum(nu);
%     % check variables that can be fixed to 0
%     indicator = lowerbound + upsilon -optval-1e-10 > 0;
%     ind = find(indicator);
%     fix0vec_final = union(fix0vec_final,ind);
%     % check variables that can be fixed to 1
%     indicator = lowerbound + nu -optval-1e-10 > 0;
%     ind = find(indicator);
%     fix1vec_final = union(fix1vec_final,ind);
% end
% 
% info.fix0vec_final = fix0vec_final;
% info.fix1vec_final = fix1vec_final;
% 
% info
