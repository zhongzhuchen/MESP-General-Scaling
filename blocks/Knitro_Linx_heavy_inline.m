%% obtain problem data
C=obj.C;
n=obj.size;
A_data=obj.A;
b_data=obj.b;
[m,~] = size(obj.A);
LB = obj.obtain_lb(s);
info = struct;

%% calling knitro
[knitro_fval,x,ininfo] = Knitro_Linx_heavy(x0,C,s,A_data,b_data,Gamma, LB);

%% assign values to info
% record important information
info.exitflag=ininfo.exitflag;
info.x=x; % optimal solution
[fval,dx,finalinfo] = obj.Linx_obj(x,s,Gamma);
info.dx=dx;
info.knitro_fval = knitro_fval;
info.fval=fval;
info.continuous_dualgap=finalinfo.dualgap;
info.dualbound=finalinfo.dualbound;
info.dual_upsilon=finalinfo.dual_upsilon;
info.dual_nu=finalinfo.dual_nu;
info.dual_pi=finalinfo.dual_pi;
info.dual_tau=finalinfo.dual_tau;
info.ub_lambda=ininfo.lambda.upper;
info.lb_lambda=ininfo.lambda.lower;

info.time=ininfo.time;
info.cputime=ininfo.cputime;
info.iterations=ininfo.output.iterations;
info.funcCount=ininfo.output.funcCount;
info.firstorderopt=ininfo.output.firstorderopt;
info.constrviolation=ininfo.output.constrviolation;
info.algorithm=ininfo.output.algorithm;

%%
fixto0list = ininfo.fixto0list;
fixto1list = ininfo.fixto1list;
info.integrality_gap=info.dualbound-obj.obtain_lb(s);
if info.integrality_gap>1e-6
    info.solved=0;
else
    info.solved=1;
end
for i=1:n
    if info.integrality_gap<info.dual_upsilon(i)-1e-10 % fix to zero
        fixto0list(end+1)=i;
    elseif info.integrality_gap<info.dual_nu(i)-1e-10 % fix to one
        fixto1list(end+1)=i;
    end
end

fixto0list = unique(fixto0list);
fixto1list = unique(fixto1list);
%% implement conflict matrices here (not the most efficient way):
%{
We store four n-by-n matrices here, where if the (i,j) element of:
1. the 1st matrix is one, then x_i and x_j cannot be one/ one simultaneously
2. the 2nd matrix is one, then x_i and x_j cannot be one/ zero simultaneously
3. the 3rd matrix is one, then x_i and x_j cannot be zero/ one simultaneously
4. the 4th matrix is one, then x_i and x_j cannot be zero/ zero simultaneously
%}
A11 = (info.integrality_gap < info.dual_upsilon + info.dual_upsilon' -1e-10);
A10 = (info.integrality_gap < info.dual_upsilon + info.dual_nu' -1e-10);
A01 = (info.integrality_gap < info.dual_nu + info.dual_upsilon' -1e-10);
A00 = (info.integrality_gap < info.dual_nu + info.dual_nu' -1e-10);

% leverage the above conflict matrices to induce more variable fixing
while true
    oldlength0 = length(fixto0list);
    oldlength1 = length(fixto1list);
    % fix more variables to 0
    if oldlength1 > 0
        A11sub_colsum = sum(A11(:, fixto1list), 2) > 0;
        fixto0list = union(fixto0list, find(A11sub_colsum));
    end
    if oldlength0 > 0
        A10sub_colsum = sum(A10(:, fixto0list), 2) > 0;
        fixto0list = union(fixto0list, find(A10sub_colsum));
    end
    % fix more variables to 0
    if oldlength1 > 0
        A01sub_colsum = sum(A01(:, fixto1list), 2) > 0;
        fixto1list = union(fixto1list, find(A01sub_colsum));
    end
    if oldlength0 > 0
        A00sub_colsum = sum(A00(:, fixto0list), 2) > 0;
        fixto1list = union(fixto1list, find(A00sub_colsum));
    end
    if oldlength0 == length(fixto0list) && oldlength1 == length(fixto1list)
        break;
    end
end

fixto0list = unique(fixto0list);
fixto1list = unique(fixto1list);

info.fix1num = length(fixto1list);
info.fixto1list = fixto1list;

info.fix0num = length(fixto0list);
info.fixto0list = fixto0list;

info.fixnum = info.fix0num + info.fix1num;


