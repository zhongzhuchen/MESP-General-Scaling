load('data63.mat');
n=length(C);
a=MESP(C,double.empty(0,n),double.empty(0,1));
for i=1:10
    for s = 40 %2:(n-1)
        scale = randi([1,3]);
        GammaInit = exp(min(max(randn(n,1),-1),1)*scale);
        [optGamma,info]=a.BFGS_DDFact_Gamma(s,GammaInit);
        if norm(optGamma-ones(n,1))>1e-5
            optGamma
        end
    end 
end
