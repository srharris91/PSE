function c=two(N);
%
% This program determines the two norm weight matrix c for
% Chebyshev polynomials. The matrix c is defined by
%
% c_{ij}= int_{-l}-{l} T_{i}(x) T_{j}(x) dx,
%
% where T_k is a Chebyshev polynomial. The above product
% satisfies
%
% c_{ij}= 1/(1-(i+j)-2)+1/(1-(i-j)-2)   for i+j even 
%                                       for i+j odd
% = 0
%
% Input
%
% N = number of modes
%
% The maximum degree M of polynomials c_i,c_j satisfies
% M = N-1
%
num=round(abs(N));
c=zeros(num,num);
for i=(0:num-1)
    for j=(0:num-1)
        if rem(i+j,2)==0,
            p=1/(1-(i+j)^2)+1/(1-(i-j)^2);
            c(i+1,j+1)=p;
        else
            c(i+1, j+1)=0;
        end;
    end;
end;
