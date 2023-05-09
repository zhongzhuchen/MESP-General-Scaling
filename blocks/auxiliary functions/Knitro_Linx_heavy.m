function [fval,x,info] = Knitro_Linx_heavy(x0,C,s,A,b,Gamma, LB)
%% obtain class properties
A_data=A;
b_data=b;
[m,~] = size(A);
scaleC=diag(Gamma)*C;
n = length(x0);
info = struct;
fixto0list = [];
fixto1list = [];
%% define nested objective function for Knitro
    function terminate = outfun(x,optimValues,state)
        if ~isinf(optimValues.fval)
            F=scaleC*diag(x)*scaleC'+eye(n)-diag(x);
            F=0.5*(F+F');
            [R,flag]=chol(F);
            Rinv=inv(R);
            K=scaleC'*Rinv;
            % calculate the derivative: 1/2*diag(C'*F^{-1}*C-F^{-1})
            dx2=sum(Rinv.*Rinv,2);
            dx=0.5*(sum(K.*K,2)-dx2)-log(Gamma);

            fN=[zeros(n,1);ones(n,1);b_data;s];
            AeqN=[-eye(n),eye(n),A_data',ones(n,1)];
            beqN=dx;
            lbN=[zeros(2*n+m,1);-inf];
            ubN=Inf(2*n+m+1,1);
            x0N=[];
            optionsN = knitro_options('algorithm',3,...  % active-set/simplex algorithm
                                     'outlev',0);       % iteration display
            [xlp, dualgap, ~, ~] = knitro_lp (fN, [], [], AeqN, beqN, lbN, ubN, x0N, [], optionsN);

            dual_upsilon = xlp(1:n);
            dual_nu = xlp((n+1):2*n);

            dualgap=dualgap+0.5*sum(dx2)-n/2+sum(x.*log(Gamma));
            integrality_gap = -optimValues.fval+ dualgap - LB;
            fixto0listN = [];
            fixto1listN = [];
            for i=1:n
                if integrality_gap<dual_upsilon(i)-1e-10 % fix to zero
                    fixto0listN(end+1)=i;
                elseif integrality_gap<dual_nu(i)-1e-10 % fix to one
                    fixto1listN(end+1)=i;
                end
            end
            fixto0list = union(fixto0list, fixto0listN);
            fixto1list = union(fixto1list, fixto1listN);
            fixto0list = unique(fixto0list);
            fixto1list = unique(fixto1list);
        end
        terminate = false;
    end
%% calling knitro
logGamma = log(Gamma);
obj_fn =  @(x) Linx_obj_Knitro_prescale(x,scaleC,logGamma);
lb=zeros(n,1);
ub=ones(n,1);
Aeq=ones(1,n);
beq=s;
A=A_data;
b=b_data;

if sum(abs(Aeq*x0-beq))>n*1e-8
    error('The initial point x0 is not feasible.')
end

TStart=tic;
tStart=cputime;

extendedFeatures.HessFcn = @(x,lambda) Linx_hessfun(x,lambda,scaleC);
extendedFeatures.OutputFcn = @ outfun;
options = knitro_options('algorithm', 0, 'convex', 1, 'derivcheck', 0, 'outlev', 0 , 'gradopt', 1, ...
                         'hessopt', 1, 'maxit', 1000, 'xtol', 1e-15, 'derivcheck_tol',1e-5,...
                         'feastol', 1e-10, 'opttol', 1e-10, 'bar_feasible',1,...
                         'bar_maxcrossit', 10);
[x,knitro_fval,exitflag,output,lambda,~] = knitro_nlp(obj_fn,x0,A,b,Aeq,beq,lb,ub,[],extendedFeatures,options);
time=toc(TStart);
tEnd=cputime-tStart;
% ========================================
fval=-knitro_fval;
info.exitflag=exitflag;
info.output=output;
info.lambda=lambda;
info.fval=fval;
info.time=time;
info.cputime = tEnd;

info.fixto0list = fixto0list;
info.fixto1list = fixto1list;
end
