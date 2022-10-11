load('data63.mat');
Corder = length(C);
Coriginal = C;
n = Corder;
% A = double.empty(0,n);
% b = double.empty(0,1);
% MESPInstance = MESP(C,A,b);
fix63 = struct;
for outs= 2:(Corder-1)
    fprintf("s = %d: \n", outs)
    allfix0=[];
    allfix1=[];
    % storing the undetermined indices
    leftind=(1:Corder)';
    s=outs;
    iteration = 0;
    C=Coriginal;
    nearzero_intgap = 0;
    while true
        n = length(C);
        x0=s/n*ones(n,1);
        A = double.empty(0,n);
        b = double.empty(0,1);
        MESPInstance = MESP(C,A,b);
        [optgamma,~]= MESPInstance.Newton_Linx_gamma(s);
        [fval,x,info_Linx] = MESPInstance.Knitro_Linx(x0,s,sqrt(optgamma)*ones(n,1));
        [fval,x,info_Fact] = MESPInstance.Knitro_DDFact(x0,s,ones(n,1));
        [fval,x,info_cFact] = MESPInstance.Knitro_DDFact_comp(x0,s,ones(n,1));
        if isempty(info_Linx.fixto0list) && isempty(info_Linx.fixto1list) && ...
                isempty(info_Fact.fixto0list) && isempty(info_Fact.fixto1list) &&...
                isempty(info_cFact.fixto0list) && isempty(info_cFact.fixto1list)
            fprintf("No more variables can be fixed.\n");
            break
        else
            fprintf("Before fixing, Linx intgap is %f, DDFact intgap is %f, DDFact comp intgap is %f. \n", ...
                info_Linx.integrality_gap, info_Fact.integrality_gap, info_cFact.integrality_gap);
            if info_Linx.integrality_gap < 1e-6 || info_Fact.integrality_gap < 1e-6 || info_cFact.integrality_gap < 1e-6 
                warning("Intergrality gap too small.\n");
                nearzero_intgap = 1;
            end
            fixto0list=union(info_Linx.fixto0list, info_Fact.fixto0list);
            fixto0list=union(fixto0list, info_cFact.fixto0list);
            fixto1list=union(info_Linx.fixto1list, info_Fact.fixto1list);
            fixto1list=union(fixto1list, info_cFact.fixto1list);
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
        fprintf("After %d times fixing, the remaing C is of order %d and rank %d, s is %d, outs is %d.\n", ...
            iteration,length(C),rank(C),s, outs);
        if isempty(C)
            fprintf("zero order C");
            break;
        end
        if s==0
            allfix0 = union(allfix0, leftind);
            fprintf("We have found the optimal solution to MESP.\n")
            break;
        end
        if length(C)==s
            allfix1 = union(allfix1, leftind);
            fprintf("We have found the optimal solution to MESP.\n")
            break;
        end
        fprintf("Number of fix to 0 variables: %d\n", length(allfix0))
        fprintf("Number of fix to 1 variables: %d\n", length(allfix1))
    end
    fprintf("----------------------------------\n");
    allfix0 = sort(allfix0);
    allfix1 = sort(allfix1);
    fix63(outs).fixto0list = allfix0;
    fix63(outs).fixto1list = allfix1;
    fix63(outs).nearzero_intgap = nearzero_intgap;
end







