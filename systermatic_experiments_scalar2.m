clear all;
load('data63.mat');
n=length(C);
folder = 'old_Results';
baseFileName = strcat('data',int2str(n),'_mix_nonuniform_lincon','.xlsx');
fullFileNameexcel = fullfile(folder, baseFileName);
if exist(fullFileNameexcel, 'file')==2
  delete(fullFileNameexcel);
end

% idx=[15
%     30
%     31
%     33
%     37
%     41
%     44
%     46
%     53
%     56
%     63];
% 
% A2=zeros(5,n);
% A2(1,idx(1:2))=1;
% A2(2,idx(3:4))=1;
% A2(3,idx(5:6))=1;
% A2(4,idx(7:8))=1;
% A2(5,idx(9:10))=1;
% 
% b2=ones(5,1);

A2 = double.empty(0,n);
b2 = double.empty(0,1);

% C=Prob.C;
% C_comp=Prob.C_comp;
% n = Prob.size;
% A_data=Prob.A;
% b_data=Prob.b;
% F=Prob.F;
% Fsquare=Prob.Fsquare;
% F_comp=Prob.F_comp;
% Fsquare_comp=Prob.Fsquare_comp;
% ldetC=Prob.ldetC;

totaltimes=1;
% A=randi([1,5],5,n);
load('random_lincon.mat');
for repeattime = 1:totaltimes
    exceloutput=[];
    for s=5:n-5
        sprintf("In the current iteration, n=%d, s=%d, repeattime=%d", n, s, repeattime)
        %% initialize the MESP instance for each (m,s) pair
        
        [xind, heurval] = heur(C,s,A2,b2);
        x=zeros(n,1);
        x(xind)=1;
        b=A*x-ones(5,1);

        Prob = MESP(C,A,b);

        lb = Prob.obtain_lb(s);
    
        x0=s/n*ones(n,1);
        %% Linx data
        [gamma,info]= Prob.BFGS_Linx_gamma(s);
        [~,~,info_Linx] = Prob.Knitro_Linx(x0,s,sqrt(gamma)*ones(n,1));
       
        %% DDFact data
        [~,~,info_DDFact] = Prob.Knitro_DDFact(x0,s,ones(n,1));
    
        %% DDFact comp data
        [~,~,info_DDFact_comp] = Prob.Knitro_DDFact_comp(x0,s,ones(n,1));
        
%         info_mix1.dualbound=nan;
%         info_mix1.integrality_gap=nan;
%         info_mix1.continuous_dualgap=nan;
%         info_mix1.fixnum_to0=nan;
%         info_mix1.fixnum_to1=nan;
%         info_mix1.fixnum=nan;
%         info_mix1.alpha=nan;
%         info_mix1.time=nan;
%         info_mix1.cputime=nan;
% 
%         info_mix2.dualbound=nan;
%         info_mix2.integrality_gap=nan;
%         info_mix2.continuous_dualgap=nan;
%         info_mix2.fixnum_to0=nan;
%         info_mix2.fixnum_to1=nan;
%         info_mix2.fixnum=nan;
%         info_mix2.alpha=nan;
%         info_mix2.time=nan;
%         info_mix2.cputime=nan;
% 
%         info_mix3.dualbound=nan;
%         info_mix3.integrality_gap=nan;
%         info_mix3.continuous_dualgap=nan;
%         info_mix3.fixnum_to0=nan;
%         info_mix3.fixnum_to1=nan;
%         info_mix3.fixnum=nan;
%         info_mix3.alpha=nan;
%         info_mix3.time=nan;
%         info_mix3.cputime=nan;
        if repeattime == 1
            %% mix DDFact and/or comp
            [~,~,info_mix1] = Prob.mix_DDFact_DDFact_comp(x0,s,ones(n,1),ones(n,1));
            %% mix DDFact Linx
            [optgamma, ~] = Prob.mix_BFGS_Gamma(s, "DDFact_Linx", gamma);
            [~,~,info_mix2] = Prob.mix_DDFact_Linx(x0,s,ones(n,1),sqrt(optgamma)*ones(n,1));
        
            %% mix DDFact comp Linx
            [optgamma, ~] = Prob.mix_BFGS_Gamma(s, "DDFact_comp_Linx", gamma);
            [~,~,info_mix3] = Prob.mix_DDFact_comp_Linx(x0,s,ones(n,1),sqrt(optgamma)*ones(n,1));
        end

        exceloutput(end+1,:)=[n, s, lb,...
            info_Linx.exitflag, info_DDFact.exitflag, info_DDFact_comp.exitflag,...
            info_Linx.dualbound, info_DDFact.dualbound, info_DDFact_comp.dualbound,...
            info_mix1.dualbound, info_mix2.dualbound, info_mix3.dualbound,...
            info_Linx.integrality_gap, info_DDFact.integrality_gap, info_DDFact_comp.integrality_gap,...
            info_mix1.integrality_gap, info_mix2.integrality_gap, info_mix3.integrality_gap,...
            info_Linx.continuous_dualgap, info_DDFact.continuous_dualgap, info_DDFact_comp.continuous_dualgap,...
            info_mix1.continuous_dualgap, info_mix2.continuous_dualgap, info_mix3.continuous_dualgap,...
            info_Linx.fixnum_to0, info_DDFact.fixnum_to0, info_DDFact_comp.fixnum_to0,...
            info_mix1.fixnum_to0, info_mix2.fixnum_to0, info_mix3.fixnum_to0,...
            info_Linx.fixnum_to1, info_DDFact.fixnum_to1, info_DDFact_comp.fixnum_to1,...
            info_mix1.fixnum_to1, info_mix2.fixnum_to1, info_mix3.fixnum_to1,...
            info_Linx.fixnum, info_DDFact.fixnum, info_DDFact_comp.fixnum,...
            info_mix1.fixnum, info_mix2.fixnum, info_mix3.fixnum,...
            info_mix1.alpha, info_mix2.alpha, info_mix3.alpha,...
            info_Linx.time, info_DDFact.time, info_DDFact_comp.time,...
            info_mix1.time, info_mix2.time, info_mix3.time,...
            info_Linx.cputime, info_DDFact.cputime, info_DDFact_comp.cputime,...
            info_mix1.cputime, info_mix2.cputime, info_mix3.cputime];
    end
    if repeattime==1
        exceloutputall=exceloutput;
    else
        exceloutputall=exceloutputall+exceloutput;
    end
end
exceloutputall=exceloutputall./totaltimes;
title=["n", "s","lb",...
    "info_Linx.exitflag", "info_DDFact.exitflag", "info_DDFact_comp.exitflag",...
        "info_Linx.dualbound", "info_DDFact.dualbound", "info_DDFact_comp.dualbound",...
        "info_mix1.dualbound", "info_mix2.dualbound", "info_mix3.dualbound",...
        "info_Linx.integrality_gap", "info_DDFact.integrality_gap", "info_DDFact_comp.integrality_gap",...
        "info_mix1.integrality_gap", "info_mix2.integrality_gap", "info_mix3.integrality_gap",...
        "info_Linx.continuous_dualgap", "info_DDFact.continuous_dualgap", "info_DDFact_comp.continuous_dualgap",...
        "info_mix1.continuous_dualgap", "info_mix2.continuous_dualgap", "info_mix3.continuous_dualgap",...
        "info_Linx.fixnum_to0", "info_DDFact.fixnum_to0", "info_DDFact_comp.fixnum_to0",...
        "info_mix1.fixnum_to0", "info_mix2.fixnum_to0", "info_mix3.fixnum_to0",...
        "info_Linx.fixnum_to1", "info_DDFact.fixnum_to1", "info_DDFact_comp.fixnum_to1",...
        "info_mix1.fixnum_to1", "info_mix2.fixnum_to1", "info_mix3.fixnum_to1",...
        "info_Linx.fixnum", "info_DDFact.fixnum", "info_DDFact_comp.fixnum",...
        "info_mix1.fixnum", "info_mix2.fixnum", "info_mix3.fixnum",...
        "info_mix1.alpha", "info_mix2.alpha", "info_mix3.alpha",...
        "info_Linx.time", "info_DDFact.time", "info_DDFact_comp.time",...
        "info_mix1.time", "info_mix2.time", "info_mix3.time",...
        "info_Linx.cputime", "info_DDFact.cputime", "info_DDFact_comp.cputime",...
        "info_mix1.cputime", "info_mix2.cputime", "info_mix3.cputime"];
xlswrite(fullFileNameexcel ,title,1,'A1');
xlswrite(fullFileNameexcel ,exceloutputall,1,'A2');


