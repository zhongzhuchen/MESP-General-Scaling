diary gscaling_iterfix_Gammaitr10_124;
load('data124.mat');
Corder = length(C);
Coriginal = C;
n = Corder;
% A = double.empty(0,n);
% b = double.empty(0,1);
% MESPInstance = MESP(C,A,b);
fix63 = struct;
% Linx_gammalist_2000 = struct;
% Linx_Gammalist_2000 = struct;
% Fact_Gammalist_2000 = struct;
% cFact_Gammalist_2000 = struct;
timelimit = 3600;
linxmaxitr = 10;
start = tic;
gscaling_list = [];
for outs= 2:(Corder-1)
    fprintf("Corder = %d, s = %d: \n", Corder, outs)
    allfix0=[];
    allfix1=[];
    % storing the undetermined indices
    leftind=(1:Corder)';
    s=outs;
    iteration = 0;
    C=Coriginal;
    nearzero_intgap = 0;
    nearzero_intgap_itr = NaN;
    while true
        n = length(C);
        x0=s/n*ones(n,1);
        A = double.empty(0,n);
        b = double.empty(0,1);
        MESPInstance = MESP(C,A,b);
        [optgamma,~]= MESPInstance.Newton_Linx_gamma(s, timelimit);
        [optGamma,info_Linxg2]= MESPInstance.BFGS_Linx_Gamma(s, sqrt(optgamma)*ones(n,1), timelimit, linxmaxitr);
        [fval,x,info_Linx] = MESPInstance.Knitro_Linx(x0,s,optGamma);
%         [optGamma,info_Factg2]=MESPInstance.BFGS_DDFact_Gamma(s,ones(n,1), timelimit);
        [fval,x,info_Fact] = MESPInstance.Knitro_DDFact(x0,s,ones(n,1));
%         [optGamma,info_cFactg2]=MESPInstance.BFGS_DDFact_comp_Gamma(s,ones(n,1), timelimit);
        [fval,x,info_cFact] = MESPInstance.Knitro_DDFact_comp(x0,s,ones(n,1));
        if isempty(info_Linx.fixto0list) && isempty(info_Linx.fixto1list) && ...
                isempty(info_Fact.fixto0list) && isempty(info_Fact.fixto1list) &&...
                isempty(info_cFact.fixto0list) && isempty(info_cFact.fixto1list)
            fprintf("No more variables can be fixed.\n");
            if iteration == 0
                Linx_fixto0list_init = [];
                Linx_fixto1list_init = [];
                Fact_fixto0list_init = [];
                Fact_fixto1list_init = [];
                cFact_fixto0list_init = [];
                cFact_fixto1list_init = [];
            end
            break
        else
            fprintf("Before fixing, Linx intgap is %f, DDFact intgap is %f, DDFact comp intgap is %f. \n", ...
                info_Linx.integrality_gap, info_Fact.integrality_gap, info_cFact.integrality_gap);
            if info_Linx.integrality_gap < 1e-6 || info_Fact.integrality_gap < 1e-6 || info_cFact.integrality_gap < 1e-6 
                warning("Intergrality gap too small.\n");
                nearzero_intgap_itr = iteration;
                nearzero_intgap = 1;
            end
            fixto0list=union(info_Linx.fixto0list, info_Fact.fixto0list);
            fixto0list=union(fixto0list, info_cFact.fixto0list);
            fixto1list=union(info_Linx.fixto1list, info_Fact.fixto1list);
            fixto1list=union(fixto1list, info_cFact.fixto1list);
            if iteration == 0
                Linx_fixto0list_init = info_Linx.fixto0list;
                Linx_fixto1list_init = info_Linx.fixto1list;
                Fact_fixto0list_init = info_Fact.fixto0list;
                Fact_fixto1list_init = info_Fact.fixto1list;
                cFact_fixto0list_init = info_cFact.fixto0list;
                cFact_fixto1list_init = info_cFact.fixto1list;
            end
            allfix0=union(allfix0,leftind(fixto0list));
            allfix1=union(allfix1,leftind(fixto1list));
            indr = setdiff((1:n)', fixto0list);
            indr = sort(indr);
            if isempty(fixto1list)
                oldC=C;
                C=C(indr,indr);
                C=1/2*(C+C');
                leftind=leftind(indr);
                leftind=sort(leftind);
            else
                % question: if the Schur complement keep the index order
                fixto1list = sort(fixto1list);
                s = s-length(fixto1list);
                indn = setdiff(indr,fixto1list);
                indn = sort(indn);
                leftind = leftind(indn);
                leftind=sort(leftind);
                D = C(fixto1list,fixto1list);
                B=C(indn,fixto1list);
                A=C(indn,indn);
                oldC=C;
                C=A-B*inv(D)*B';
                C=1/2*(C+C');
            end
        end
        iteration=iteration+1;
        fprintf('-');
        fprintf("After %d times fixing, the remaing C is of order %d and rank %d, s is %d, outs is %d.\n", ...
            iteration,length(C),rank(C),s, outs);
        fprintf("BFGS iterations. Linx: %d.\n", ...
            info_Linxg2.iterations);
        fprintf("Time: %f.\n", toc(start));
        breakindicator = 0;
        if isempty(C)
            fprintf("zero order C.\n");
            breakindicator = 1;
        end
        if s==0
            allfix0 = union(allfix0, leftind);
            fprintf("We have found the optimal solution to MESP.\n")
            breakindicator = 1;
        end
        if length(C)==s
            allfix1 = union(allfix1, leftind);
            fprintf("We have found the optimal solution to MESP.\n")
            breakindicator = 1;
        end
        fprintf("Number of fix to 0 variables: %d\n", length(allfix0))
        fprintf("Number of fix to 1 variables: %d\n", length(allfix1))
        if breakindicator == 1
            break;
        end
    end
    fprintf("----------------------------------\n");
    allfix0 = sort(allfix0);
    allfix1 = sort(allfix1);
    gscaling_list(end+1)=min(Corder, length(allfix0)+length(allfix1));
    fix63(outs).fixto0list = allfix0;
    fix63(outs).fixto1list = allfix1;
    fix63(outs).Linx_fixto0list_init = Linx_fixto0list_init;
    fix63(outs).Linx_fixto1list_init = Linx_fixto1list_init;
    fix63(outs).Fact_fixto0list_init = Fact_fixto0list_init;
    fix63(outs).Fact_fixto1list_init = Fact_fixto1list_init;
    fix63(outs).cFact_fixto0list_init = cFact_fixto0list_init;
    fix63(outs).cFact_fixto1list_init = cFact_fixto1list_init;
    fix63(outs).nearzero_intgap = nearzero_intgap;
    fix63(outs).nearzero_intgap_itr = nearzero_intgap_itr;
    fix63(outs).iteration = iteration;
end
if Corder == 63
    fixnumsnew.gscaling_63 = gscaling_list;
elseif Corder == 90
    fixnumsnew.gscaling_90 = gscaling_list;
elseif Corder == 124
    fixnumsnew.gscaling_124 = gscaling_list;
elseif Corder == 150
    fixnumsnew.gscaling_150 = gscaling_list;
elseif Corder == 300
    fixnumsnew.gscaling_300 = gscaling_list;
end
% fix=fix63;
% filename = strcat('iterative_fixing/fix', num2str(Corder), '_LFcF_gscaling.mat');
% save(filename,'fix');
diary off;







