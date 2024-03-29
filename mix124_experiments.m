clear all:
load('data124.mat');
n=length(C);
folder = 'mix_Results';
baseFileName = strcat('data',int2str(n),'.xlsx');
fullFileNameexcel = fullfile(folder, baseFileName);
if exist(fullFileNameexcel, 'file')==2
  delete(fullFileNameexcel);
end
exceloutput=[];

A = double.empty(0,n);
b = double.empty(0,1);
for s=90:(n-1)
    s
    Prob = MESP(C,A,b);
    %% lower bound
    lb = Prob.obtain_lb(s);
    x0 = s/n*ones(n,1);
    %% Linx Data
    [optgamma,info]= Prob.BFGS_Linx_gamma(s);
    [~,~,info_Linx] = Prob.Knitro_Linx(x0,s,sqrt(optgamma)*ones(n,1));

    %% DDFact Data
    [~,~,info_DDFact] = Prob.Knitro_DDFact(x0,s,ones(n,1));

    %% DDFact comp data
    [~,~,info_DDFact_comp] = Prob.Knitro_DDFact_comp(x0,s,ones(n,1));

    %% mix all
    [gamma,info_mix]=Prob.mix_DDFact_DDFact_comp_Linx_BFGS_Gamma(s,optgamma);

    exceloutput(end+1,:)=[n,s,lb,...
        info_Linx.dualbound, info_DDFact.dualbound, info_DDFact_comp.dualbound, info_mix.dualbound,...
        info_Linx.integrality_gap, info_DDFact.integrality_gap, info_DDFact_comp.integrality_gap, info_mix.integrality_gap,...
        info_Linx.continuous_dualgap, info_DDFact.continuous_dualgap, info_DDFact_comp.continuous_dualgap, info_mix.continuous_dualgap,...
        info_Linx.fixnum_to0, info_DDFact.fixnum_to0, info_DDFact_comp.fixnum_to0, info_mix.fixnum_to0,...
        info_Linx.fixnum_to1, info_DDFact.fixnum_to1, info_DDFact_comp.fixnum_to1, info_mix.fixnum_to1,...
        info_Linx.fixnum, info_DDFact.fixnum, info_DDFact_comp.fixnum, info_mix.fixnum,...
        info_mix.alpha(1), info_mix.alpha(2), info_mix.alpha(3),...
        info_Linx.time, info_DDFact.time, info_DDFact_comp.time, info_mix.time,...
        info_Linx.cputime, info_DDFact.cputime, info_DDFact_comp.cputime, info_mix.cputime];
end
title=["n","s","lb",...
        "info_Linx.dualbound," "info.DDFact.dualbound," "info_DDFact_comp.dualbound," "info_mix.dualbound,"...
        "info_Linx.integrality_gap," "info.DDFact.integrality_gap," "info_DDFact_comp.integrality_gap," "info_mix.integrality_gap,"...
        "info_Linx.continuous_dualgap," "info.DDFact.continuous_dualgap," "info_DDFact_comp.continuous_dualgap," "info_mix.continuous_dualgap,"...
        "info_Linx.fixnum_to0," "info.DDFact.fixnum_to0," "info_DDFact_comp.fixnum_to0," "info_mix.fixnum_to0,"...
        "info_Linx.fixnum_to1," "info.DDFact.fixnum_to1," "info_DDFact_comp.fixnum_to1," "info_mix.fixnum_to1,"...
        "info_Linx.fixnum," "info.DDFact.fixnum," "info_DDFact_comp.fixnum," "info_mix.fixnum,"...
        "info_mix.alpha(1)," "info_mix.alpha(2)," "info_mix.alpha(3),"...
        "info_Linx.time," "info.DDFact.time," "info_DDFact_comp.time," "info_mix.time,"...
        "info_Linx.cputime," "info.DDFact.cputime," "info_DDFact_comp.cputime," "info_mix.cputime"];
 
xlswrite(fullFileNameexcel ,title,1,'A1');
xlswrite(fullFileNameexcel ,exceloutput,1,'A2');
