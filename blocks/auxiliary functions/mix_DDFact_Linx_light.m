function [fval,x,info] = mix_DDFact_Linx_light(x0,C,s,F,Fsquare,A_data,b_data,Gamma1,Gamma2)
n=length(C);
[m,~] = size(A_data);
info=struct;
%% solve linx and fact bound
TStart=tic;
tStart=cputime;
% obtain optimal solutionnn at search boundary a=0, b=1
[~, xa, ~] = Knitro_DDFact_light(x0,C,s,F,Fsquare,A_data,b_data,Gamma1);
[~, xb, ~] = Knitro_Linx_light(x0,C,s,A_data,b_data,Gamma2);
[fval_a1, dx_a1, info_a1] = DDFact_obj_Knitro(xa,s,F,Fsquare,Gamma1);
[fval_a2, dx_a2, info_a2] = Linx_obj_Knitro(xa,C,Gamma2);
[fval_b1, dx_b1, info_b1] = DDFact_obj_Knitro(xb,s,F,Fsquare,Gamma1);
[fval_b2, dx_b2, info_b2] = Linx_obj_Knitro(xb,C,Gamma2);
% gradient of the mixing objective function
dx_a = -dx_a1;
dx_b = -dx_b2;
% gradient of the mixing parameter alpha
dalpha_a = -fval_a2+fval_a1;
dalpha_b = -fval_b2+fval_b1;

if dalpha_a>=0
    fval = -fval_a1;
    x = xa;
    dx = dx_a;
    alpha = 0;
    info1 = info_a1;
    info2 = info_a2;
elseif dalpha_b <=0
    fval = -fval_b2;
    x = xb;
    dx = dx_b;
    alpha = 1;
    info1 = info_b1;
    info2 = info_b2;
else
    lb=zeros(n,1);
    ub=ones(n,1);
    Aeq=ones(1,n);
    beq=s;
    options = knitro_options('algorithm', 0, 'convex', 1, 'derivcheck', 0, 'outlev', 0 , 'gradopt', 1, ...
                             'hessopt', 2, 'maxit', 1000, 'xtol', 1e-15, ...
                             'feastol', 1e-10, 'opttol', 1e-10, 'bar_feasible',1,...
                             'bar_maxcrossit', 10);
    alpha=1/2;
    obj_fn = @(x) mix_DDFact_Linx_obj_Knitro(x,C,s,F,Fsquare,Gamma1,Gamma2,alpha);
    try
    [x,knitro_fval,exitflag,output,lambda,~] = knitro_nlp(obj_fn,x0,A_data,b_data,Aeq,beq,lb,ub,[],[],options);
    catch
        1+1
    end
    [fval1, dx1, info1] = DDFact_obj_Knitro(x,s,F,Fsquare,Gamma1);
    [fval2, dx2, info2] = Linx_obj_Knitro(x,C,Gamma2);
    fval=-knitro_fval;
    dx=-(1-alpha)*dx1-alpha*dx2;
    dalpha=-fval2+fval1;
    % use the proximal gradient method
    beta=0.5; % shrink parameter for t
    Tol=1e-8;
    Numiterations=100;
    k=0;
    while true
        k=k+1;
        t=1;
        while true
            if t<1e-6
                break;
            end
            Gt=1/t*(alpha-max(min(alpha-t*dalpha,1),0));
            newalpha=alpha-t*Gt;
            obj_fn = @(x) mix_DDFact_Linx_obj_Knitro(x,C,s,F,Fsquare,Gamma1,Gamma2,newalpha);
            % warm start
            try
            [x,knitro_fval,exitflag,output,lambda,~] = knitro_nlp(obj_fn,x0,A_data,b_data,Aeq,beq,lb,ub,[],[],options);
            catch
                pause;
            end
            [fval1, dx1, info1] = DDFact_obj_Knitro(x,s,F,Fsquare,Gamma1);
            [fval2, dx2, info2] = Linx_obj_Knitro(x,C,Gamma2);
            newfval=-knitro_fval;
            newdalpha=-fval2+fval1;
            if newfval<=fval-t*dalpha'*Gt+t/2*(Gt'*Gt)
                break
            end
            t=beta*t;
        end
        if k>Numiterations || (Gt'*Gt)<Tol ||fval-newfval<Tol
            alpha=newalpha;
            dalpha=newdalpha;
            fval=newfval;
            dx=-(1-alpha)*dx1-alpha*dx2;
            break
        end
        alpha=newalpha;
        dalpha=newdalpha;
        fval=newfval;
        dx=-(1-alpha)*dx1-alpha*dx2;
    end
end
time=toc(TStart);
tEnd=cputime-tStart;
%% assign values to info
% record important information
info.x=x; % optimal solution
info.dx=dx;
info.fval=fval;
info.alpha = alpha;
info.time=time;
info.cputime=tEnd;
end

