%% obtain problem data
info=struct;
C=obj.C;
n=obj.size;
A_data=obj.A;
b_data=obj.b;
[m,~] = size(obj.A);
F_comp=obj.F_comp;
Fsquare_comp = obj.Fsquare_comp;
ldetC = obj.ldetC;

%% calling knitro to solve the DDFact relaxation problem
[knitro_fval,x,ininfo] = Knitro_DDFact_comp_light(x0,C,s,F_comp,Fsquare_comp,ldetC,A_data,b_data,Gamma);

%% assign values to info
% record important information
info.exitflag=ininfo.exitflag;
info.x=x; % optimal solution
[fval,dx,finalinfo] = obj.DDFact_comp_obj(x,s,Gamma);
info.dx=dx;
info.knitro_fval = knitro_fval;
info.fval=fval;
info.continuous_dualgap=finalinfo.dualgap;
info.dualbound=finalinfo.dualbound;
% Please note that the dual variables here are associate with y=e-x
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

%% fixing variables
info.fixnum=0;
info.fixnum_to0=0;
info.fixto0list=[];
info.fixnum_to1=0;
info.fixto1list=[];

% if the integrality gap between the upper bound and the lower bound is
% large, the fixing power might be hidden due the weak lower bound
info.integrality_gap=info.dualbound-obj.obtain_lb(s);
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

