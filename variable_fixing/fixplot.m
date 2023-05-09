load('data63fix.mat');
load('data63solution.mat');
n = length(C);
fixM0 = zeros(n);
fixM1 = zeros(n);
sol0 = zeros(n);
sol1 = zeros(n);

for s =2:(n-1)
    fieldname0 = strcat('g_0_', num2str(s));
    fieldname1 = strcat('g_1_', num2str(s));
    vec0 = zeros(1,n);
    vec0(info_maxfix.(fieldname0)) = 1;
    fixM0(s, :) = vec0;
    fieldname = strcat('x', num2str(s));
    sol0(s,:) = 1-info_x.(fieldname);

    vec1 = zeros(1,n);
    vec1(info_maxfix.(fieldname1)) = 1;
    fixM1(s, :) = vec1;
    sol1(s,:) = 1-info_x.(fieldname);
end
% figure;
% imagesc(fixM0);
% axis on;
% xlabel("element");
% ylabel("s");
% figure;
% plot(sol0);zhongz
figure;
imagesc(fixM1);
axis on;
xlabel("element");
ylabel("s");
figure;
scatter(1:n, sol1');

% load('data63fix.mat');
% n = length(C);
% info_x = struct;
% for s=2:62
%     [optgamma,info]= MESPEX1.Newton_Linx_gamma(s, 300);
%     [optGamma,info]= MESPEX1.BFGS_Linx_Gamma(s, sqrt(optgamma)*ones(n,1), 300);
%     [fval,x,info] = MESPEX1.Knitro_Linx(s/n*ones(n,1),s,optGamma);
%     fieldname = strcat('x', num2str(s));
%     info_x.(fieldname) = info.x;
% end