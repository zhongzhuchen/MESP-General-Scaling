% load('data63.mat');
% n=length(C);
% a=MESP(C,double.empty(0,n),double.empty(0,1));
% for i=1:10
%     for s = 40 %2:(n-1)
%         scale = randi([1,3]);
%         GammaInit = exp(min(max(randn(n,1),-1),1)*scale);
%         [optGamma,info]=a.BFGS_DDFact_Gamma(s,GammaInit);
%         if norm(optGamma-ones(n,1))>1e-5
%             optGamma
%         end
%     end 
% end

% s = 20: 
% Before fixing, Linx intgap is 1.100141, DDFact intgap is 0.507719, DDFact comp intgap is 4.139407. 
% After 1 times fixing, the remaing C is of order 66 and rank 66, s is 20, outs is 20.
% Number of fix to 0 variables: 58
% Number of fix to 1 variables: 0
% No more variables can be fixed.
% ----------------------------------
% s = 40: 
% No more variables can be fixed.
% ----------------------------------
% s = 60: 
% No more variables can be fixed.
% ----------------------------------
% s = 80: 
% No more variables can be fixed.
% ----------------------------------
% s = 100: 
% No more variables can be fixed.
% ----------------------------------
% iterative_fixing_singular
% s = 20: 
% Before fixing, Linx intgap is 1.119965, DDFact intgap is 0.101430. 
% After 1 times fixing, the remaing C is of order 28 and rank 28, s is 20, outs is 20.
% Number of fix to 0 variables: 1972
% Number of fix to 1 variables: 0
% Before fixing, Linx intgap is 0.016507, DDFact intgap is 0.101430. 
% After 2 times fixing, the remaing C is of order 10 and rank 10, s is 2, outs is 20.
% Number of fix to 0 variables: 1972
% Number of fix to 1 variables: 18
% Before fixing, Linx intgap is 0.000000, DDFact intgap is -0.000000. 
% Warning: Intergrality gap too small.\n 
% > In iterative_fixing_singular (line 53) 
% After 3 times fixing, the remaing C is of order 1 and rank 1, s is 0, outs is 20.
% We have found the optimal solution to MESP.
% ----------------------------------
% s = 40: 
% Before fixing, Linx intgap is 2.753006, DDFact intgap is 0.235859. 
% After 1 times fixing, the remaing C is of order 58 and rank 58, s is 40, outs is 40.
% Number of fix to 0 variables: 1942
% Number of fix to 1 variables: 0
% Before fixing, Linx intgap is 0.056662, DDFact intgap is 0.235859. 
% After 2 times fixing, the remaing C is of order 23 and rank 23, s is 6, outs is 40.
% Number of fix to 0 variables: 1943
% Number of fix to 1 variables: 34
% Before fixing, Linx intgap is 0.018422, DDFact intgap is 0.039717. 
% After 3 times fixing, the remaing C is of order 6 and rank 6, s is 3, outs is 40.
% Number of fix to 0 variables: 1957
% Number of fix to 1 variables: 37
% Before fixing, Linx intgap is 0.001338, DDFact intgap is 0.018342. 
% After 4 times fixing, the remaing C is of order 4 and rank 4, s is 2, outs is 40.
% Number of fix to 0 variables: 1958
% Number of fix to 1 variables: 38
% Before fixing, Linx intgap is 0.000287, DDFact intgap is 0.003891. 
% After 5 times fixing, the remaing C is of order 2 and rank 2, s is 1, outs is 40.
% Number of fix to 0 variables: 1959
% Number of fix to 1 variables: 39
% Before fixing, Linx intgap is 0.000000, DDFact intgap is 0.000000. 
% Warning: Intergrality gap too small.\n 
% > In iterative_fixing_singular (line 53) 
% After 6 times fixing, the remaing C is of order 1 and rank 1, s is 1, outs is 40.
% We have found the optimal solution to MESP.
% ----------------------------------
% s = 60: 
% Before fixing, Linx intgap is 5.806206, DDFact intgap is 0.303972. 
% After 1 times fixing, the remaing C is of order 110 and rank 107, s is 60, outs is 60.
% Number of fix to 0 variables: 1890
% Number of fix to 1 variables: 0
% Before fixing, Linx intgap is 0.130397, DDFact intgap is 0.303972. 
% After 2 times fixing, the remaing C is of order 67 and rank 64, s is 17, outs is 60.
% Number of fix to 0 variables: 1890
% Number of fix to 1 variables: 43
% Before fixing, Linx intgap is 0.065952, DDFact intgap is 0.095332. 
% After 3 times fixing, the remaing C is of order 18 and rank 18, s is 11, outs is 60.
% Number of fix to 0 variables: 1933
% Number of fix to 1 variables: 49
% Before fixing, Linx intgap is 0.000000, DDFact intgap is 0.024896. 
% Warning: Intergrality gap too small.\n 
% > In iterative_fixing_singular (line 53) 
% After 4 times fixing, the remaing C is of order 1 and rank 1, s is 0, outs is 60.
% We have found the optimal solution to MESP.
% ----------------------------------
% s = 80: 
% Before fixing, Linx intgap is 8.986346, DDFact intgap is 0.612019. 
% After 1 times fixing, the remaing C is of order 442 and rank 286, s is 80, outs is 80.
% Number of fix to 0 variables: 1558
% Number of fix to 1 variables: 0
% No more variables can be fixed.
% ----------------------------------
% s = 100: 
% No more variables can be fixed.
