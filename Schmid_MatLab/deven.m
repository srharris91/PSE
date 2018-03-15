function d1=deven(N);
%
% Compute matrix which converts Chebyshev coefficients of a
% polynomial to coefficients of derivative
%
% Reference:
% Gottlieb and Orszag, Numerical Analysis of Spectral
% Methods: Theory and Applications, SIAM, Philadelphia,
% 1977.
%
% d1 derivative matrix (N,N)
% N number of coefficients
%
num=round(abs(N));
d1=zeros(num,num);
for i=0: (num-1)
    for j=(i+1):2:(num-1)
        d1(i+1,j+1)=2*j;
    end;
end;
d1(1,:)=d1(1,:)/2;

