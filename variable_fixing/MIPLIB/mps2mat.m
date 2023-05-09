function mps2mat(problem_name, optval)
% a function to convert mps file into .mat
problem_name_mps = strcat(problem_name, '.mps');
problem = mpsread(problem_name_mps);
f = problem.f;
A = problem.Aineq;
b = problem.bineq;
Aeq = problem.Aeq;
beq = problem.beq;
lb = problem.lb;
ub = problem.ub;
integrality = problem.intcon;
matname = strcat(problem_name, '.mat');
save(matname, 'f', 'A', 'b', 'Aeq', 'beq', 'lb', 'ub', "integrality", 'optval');
sprintf('saved');
end

