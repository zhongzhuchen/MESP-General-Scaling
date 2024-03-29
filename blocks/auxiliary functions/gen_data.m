function [F,Fsquare,ldetC] = gen_data(C,comp)
%generate the decomposition C= FF' for the use of calculating DDFact bound
%{
Input:
    C         - data matrix
    comp      - equals to 1 means we decompose C; equals to 0 means we
                decompose inv(C)

Output:
    F         - a factor of C such that C=F*F', in the dimension n-by-d
    Fsquare   - a 3d array where Fsquare(:,:,i) represents the F(i,:)'*F(i,:)
    ldetC     - logarithm determinant of C
%}
n=length(C);
ldetC=log(det(C));
if comp==1
    C=inv(C);
end
d=rank(C);
[U,D]=eig(C);
s=diag(D);
[s,I] = sort(s, 'descend');

% %% force eigenvalues corresponding to position > d = rank(C) to be zero
% s(d+1:end) = 0;

U = U(:, I);
F=diag(sqrt(s))*U';

thresh = 0;
if d + thresh <= n
    F=F(1:(d+thresh),1:n);
    F=F';
    Fsquare=zeros(d+thresh,d+thresh,n);
end

for i=1:n
    Fsquare(:,:,i)=F(i,:)'*F(i,:);
end
end

