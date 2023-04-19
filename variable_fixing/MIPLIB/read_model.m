% problem = mpsread('presolved/academictimetablesmall_presolved.mps'); optval = 0;
problem = mpsread('mps/academictimetablesmall.mps'); optval = 0;
% problem = mpsread('mps/bnatt400.mps'); optval = 1;
% problem = mpsread('presolved/bnatt400_presolved.mps'); optval = 1;
% problem = mpsread('mps/neos-3555904-turama.mps'); optval = -34.7;
% [info] = Knitro_LP_heavy(problem, optval);
% sprintf("===============================")
% [Dinfo] = Knitro_DLP_heavy(problem, optval);
% sprintf("===============================")
% [Ginfo] = Gurobi_fix(problem, optval);
% [fix0vec, fix1vec, info] = maxfix(problem, optval);
