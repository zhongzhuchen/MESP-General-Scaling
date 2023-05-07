function [fval,x,info]=Knitro_DDFact_heavy(x0,C,s,F,Fsquare,A,b,Gamma, LB, param)
%% obtain problem data 
A_data=A;
b_data=b;
[m,~] = size(A);
scaleC=diag(Gamma)*C;
n = length(x0);
d=length(Fsquare(:,:,1));
info = struct;
info.boundarycount1 = 0;
info.boundarycount2 = 0;
fixto0list = [];
fixto1list = [];
%% define nested objective function for Knitro
    function terminate = outfun(x,optimValues,state)
%         if ~isinf(optimValues.fval)
%             [fval,dx,~] = DDFact_obj_Knitro(x,s,F,Fsquare,Gamma);
%             dx=-dx;
%             fval=-fval;
%             fN=[zeros(n,1);ones(n,1);b_data;s];
%             AeqN=[-eye(n),eye(n),A_data',ones(n,1)];
%             beqN=dx;
%             lbN=[zeros(2*n+m,1);-inf];
%             ubN=Inf(2*n+m+1,1);
%             x0N=[];
%             optionsN = knitro_options('algorithm',3,...  % active-set/simplex algorithm
%                                      'outlev',0);       % iteration display
%             [xlp, dualgap, ~, ~] = knitro_lp (fN, [], [], AeqN, beqN, lbN, ubN, x0N, [], optionsN);
% 
%             dual_upsilon = xlp(1:n);
%             dual_nu = xlp((n+1):2*n);
% 
%             dualgap=dualgap-s+sum(x.*log(Gamma));
%             integrality_gap = -optimValues.fval+ dualgap - LB;
%             fixto0listN = [];
%             fixto1listN = [];
%             for i=1:n
%                 if integrality_gap<dual_upsilon(i)-1e-10 % fix to zero
%                     fixto0listN(end+1)=i;
%                 elseif integrality_gap<dual_nu(i)-1e-10 % fix to one
%                     fixto1listN(end+1)=i;
%                 end
%             end
%             fixto0list = union(fixto0list, fixto0listN);
%             fixto1list = union(fixto1list, fixto1listN);
%             fixto0list = unique(fixto0list);
%             fixto1list = unique(fixto1list);
%         end
        % check if x locates on the boundary of the function domain
        X=zeros(d);
        for i=1:n
            if x(i)==0
            else
                X=X+x(i)*Fsquare(:,:,i);
            end
        end
%         [~,D] = eig(X);
%         D = sort(diag(D),'descend');
%         if D(end) < 1e-12 
%             info.boundarycount1 = info.boundarycount1 + 1;
%         end
        if rank(X) < d
            info.boundarycount1 = info.boundarycount1 + 1;
        end
        if min(x) < 1e-14
            info.boundarycount2 = info.boundarycount2 + 1;
        end
%         sprintf("the minimum eigenvalue: %0.15e, minimum value of x: %e", D(end), min(x));
        terminate = false;
    end
%% calling knitro
obj_fn =  @(x) DDFact_obj_Knitro(x,s,F,Fsquare,Gamma);
lb=zeros(n,1);
ub=ones(n,1);
Aeq=ones(1,n);
beq=s;
A=A_data;
b=b_data;

if sum(abs(Aeq*x0-beq))>n*1e-10
    error('The initial point x0 is not feasible.')
end

extendedFeatures.OutputFcn = @ outfun;
options = knitro_options('algorithm', param.algorithm, 'convex', 1, 'derivcheck', 0, 'outlev', 0 , 'gradopt', 1, ...
                         'hessopt', 2, 'maxit', 1000, 'xtol', 1e-15, ...
                         'feastol', 1e-10, 'opttol', 1e-10, 'bar_feasible',1,...
                         'bar_maxcrossit', 10);
TStart=tic;
tStart=cputime;
[x,knitro_fval,exitflag,output,lambda,~] = knitro_nlp(obj_fn,x0,A,b,Aeq,beq,lb,ub,[],extendedFeatures,options);  
time=toc(TStart);
tEnd=cputime-tStart;

fval=-knitro_fval;
info.exitflag=exitflag;
info.output=output;
info.lambda=lambda;
info.fval=fval;
info.time=time;
info.cputime = tEnd;
end
