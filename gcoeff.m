function PHIm = gcoeff(lam,N,K)
% lam=[ 0 0 0 0 0 1 2 3 5 7 ];
% N=2;
% K=2;
% ;
% N=5;
% K=3;


M=length(lam);
c=M/N;

B='0';
for i=1:M
    B=[B '+' num2str(lam(i)) '/(' num2str(lam(i)) '-x)' ];
end
B=[B '-' num2str(N)];
syms x;
A=solve(B);
miu0=(sort(double(A)));

n=length(miu0);

miu=zeros(1,M);

miu(M-n+1:M)=miu0;
PHIm=zeros(1,M);
for m=1:M
    if m<=M-K
        k=M-K+1:M;
        PHIm(m)=1+sum(lam(k)./(lam(m)-lam(k))-miu(k)./(lam(m)-miu(k)));
    else 
        k=1:M-K;
        PHIm(m)=-sum(lam(k)./(lam(m)-lam(k))-miu(k)./(lam(m)-miu(k)));
    end
end
