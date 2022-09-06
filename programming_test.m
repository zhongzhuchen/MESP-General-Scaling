load('data63.mat');
load('random_lincon.mat');
s=31;
n=length(C);

% constructor check
MESPInstance = MESP(C,A,b);
fprintf('<MESP> constructer method passed.\n');

F = MESPInstance.F;
Fsquare = MESPInstance.Fsquare;

F_comp = MESPInstance.F_comp;
Fsquare_comp = MESPInstance.Fsquare_comp;

if norm(F*F'-C)>1e-8 & norm(F_comp*F_comp'-inv(C))>1e-8
    error('<gen_data> decomposition method failed.\n');
else
    fprintf('<gen_data> decomposition method passed.\n');
end

% % heuristic solution + lower bound check
% [lb,info] = MESPInstance.obtain_lb(s);
% if not(isnan(info.x))
%     if length(info.x)>n | length(info.x)<0
%         error('Heuristic method <heur> failed.\n');
%     else
%         fprintf('Heuristic method <heur> passed.\n');
%     end
% end
% 
% % factorization bound check
% x0 = s/n*ones(n,1);
% Gamma = ones(n,1);
% [fval,x,info] = MESPInstance.Knitro_DDFact(x0,s,Gamma);
% if abs(sum(x)-s)>1e-10 | min(x)<-1e-10 | max(x)>1+1e-10
%     error('<Knitro_DDFact> method does not give a final feasible solution.\n')
% elseif info.continuous_dualgap>1e-8
%     error('<Knitro_DDFact> method does not give an final optimal solution.\n')
% end
% 
% Gamma = rand(n,1)+0.5;
% [fval,x,info] = MESPInstance.Knitro_DDFact(x0,s,Gamma);
% if abs(sum(x)-s)>1e-10 | min(x)<-1e-10 | max(x)>1+1e-10
%     error('<Knitro_DDFact> method does not give a final feasible solution.\n')
% elseif info.continuous_dualgap>1e-8
%     error('<Knitro_DDFact> method does not give an final optimal solution.\n')
% else
%     fprintf('<Knitro_DDFact> method passed.\n')
% end
% 
% % general scaling vector for factorization optimization check
% GammaInit = ones(n,1);
% [optGamma,info] = MESPInstance.BFGS_DDFact_Gamma(s,GammaInit);
% if info.absres>1e-6
%     error('<BFGS_DDFact_Gamma> method does not give a local optimal solution.\n')
% else
%     fprintf('<BFGS_DDFact_Gamma> method passed.\n')
% end
% 
% % complementary factorization bound check
% x0 = s/n*ones(n,1);
% Gamma = ones(n,1);
% [fval,x,info] = MESPInstance.Knitro_DDFact_comp(x0,s,Gamma);
% if abs(sum(x)-s)>1e-10 | min(x)<-1e-10 | max(x)>1+1e-10
%     error('<Knitro_DDFact_comp> method does not give a final feasible solution.\n')
% elseif info.continuous_dualgap>1e-8
%     error('<Knitro_DDFact_comp> method does not give an final optimal solution.\n')
% end
% 
% Gamma = rand(n,1)+0.5;
% [fval,x,info] = MESPInstance.Knitro_DDFact_comp(x0,s,Gamma);
% if abs(sum(x)-s)>1e-10 | min(x)<-1e-10 | max(x)>1+1e-10
%     error('<Knitro_DDFact_comp> method does not give a final feasible solution.\n')
% elseif info.continuous_dualgap>1e-8
%     error('<Knitro_DDFact_comp> method does not give an final optimal solution.\n')
% else
%     fprintf('<Knitro_DDFact_comp> method passed.\n')
% end
% 
% % general scaling vector for complementary factorization optimization check
% GammaInit = ones(n,1);
% [optGamma,info] = MESPInstance.BFGS_DDFact_comp_Gamma(s,GammaInit);
% if info.absres>1e-6
%     error('<BFGS_DDFact_comp_Gamma> method does not give a local optimal solution.\n')
% else
%     fprintf('<BFGS_DDFact_comp_Gamma> method passed.\n')
% end

% Linx bound check
x0 = s/n*ones(n,1);
Gamma = ones(n,1);
[fval,x,info] = MESPInstance.Knitro_Linx(x0,s,Gamma);
if abs(sum(x)-s)>1e-10 | min(x)<-1e-10 | max(x)>1+1e-10
    error('<Knitro_Linx> method does not give a final feasible solution.\n')
elseif info.continuous_dualgap>1e-8
    error('<Knitro_Linx> method does not give an final optimal solution.\n')
end

Gamma = rand(n,1)+0.5;
[fval,x,info] = MESPInstance.Knitro_Linx(x0,s,Gamma);
if abs(sum(x)-s)>1e-10 | min(x)<-1e-10 | max(x)>1+1e-10
    error('<Knitro_Linx> method does not give a final feasible solution.\n')
elseif info.continuous_dualgap>1e-8
    error('<Knitro_Linx> method does not give an final optimal solution.\n')
else
    fprintf('<Knitro_Linx> method passed.\n')
end