%% obtain class properties
A_data=obj.A;
b_data=obj.b;
[m,~] = size(obj.A);
info=struct;
scaleC=diag(Gamma)*obj.C;
n = obj.size;
F_comp = obj.F_comp;
Fsquare_comp = obj.Fsquare_comp;
ldetC = obj.ldetC;
%% calling knitro
obj_fn =  @(x) DDFact_comp_obj_Knitro(x,s,F_comp,Fsquare_comp,ldetC,Gamma);
lb=zeros(n,1);
ub=ones(n,1);
Aeq=ones(1,n);
beq=s;
A=A_data;
b=b_data;

if sum(abs(Aeq*x0-beq))>1e-10
    error('The initial point x0 is not feasible.')
end

options = knitro_options('algorithm', 3, 'convex', 1, 'derivcheck', 0, 'outlev', 0 , 'gradopt', 1, ...
                         'hessopt', 2, 'maxit', 1000, 'xtol', 1e-15,...
                         'feastol', 1e-10, 'opttol', 1e-10, 'bar_feasible',1,...
                         'bar_maxcrossit', 10);
TStart=tic;
tStart=cputime;
[x,knitro_fval,exitflag,output,lambda,~] = knitro_nlp(obj_fn,x0,A,b,Aeq,beq,lb,ub,[],[],options);  
time=toc(TStart);
tEnd=cputime-tStart;

%% assign values to info
% record important information
info.exitflag=exitflag;
info.x=x; % optimal solution
info.knitro_fval = -knitro_fval;
[fval,dx,finalinfo] = obj.DDFact_comp_obj(x,s,Gamma);
info.dx=dx;
info.fval=fval;
info.continuous_dualgap=finalinfo.dualgap;
info.dualbound=finalinfo.dualbound;
% Please note that the dual variables here are associate with y=e-x
info.dual_upsilon=finalinfo.dual_upsilon;
info.dual_nu=finalinfo.dual_nu;
info.dual_pi=finalinfo.dual_pi;
info.dual_tau=finalinfo.dual_tau;

% Please note that these variables here are associate with x
info.ub_lambda=lambda.upper;
info.lb_lambda=lambda.lower;

info.time=time;
info.cputime=tEnd;
info.iterations=output.iterations;
info.funcCount=output.funcCount;
info.firstorderopt=output.firstorderopt;
info.constrviolation=output.constrviolation;
info.algorithm=output.algorithm;

%% fixing variables
info.fixnum=0;
info.fixnum_to0=0;
info.fixto0list=[];
info.fixnum_to1=0;
info.fixto1list=[];

info.integrality_gap=info.dualbound-obj.obtain_lb(s);
if info.integrality_gap>1e-6
    info.solved=0;
else
    info.solved=1;
end
for i=1:n
    if info.integrality_gap<info.dual_upsilon(i)-1e-10
        info.fixnum=info.fixnum+1;
        info.fixnum_to0=info.fixnum_to0+1;
        info.fixto0list(end+1)=i;
    elseif info.integrality_gap<info.dual_nu(i)-1e-10
        info.fixnum=info.fixnum+1;
        info.fixnum_to1=info.fixnum_to1+1;
        info.fixto1list(end+1)=i;
    end
end

