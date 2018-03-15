function t=fmy(t1,t2,options);
%
% This function uses the built-in function 'fmin' to find
% the maximum value of a function on the interval [t1,t2].
% The function is in the file maxer.m
%
% INPUT:
% t1, t2 lower and upper bounds of interval
% options input parameters for minimization routine
%
% OUTPUT
% t = value at which function maxer(t) is minimized
f1=maxer(t1);
f2=maxer(t2);
tt=fminbnd('maxer' ,t1,t2,options);
f3=maxer(tt);
f=[f1 f2 f3];
tm=[t1 t2 tt];
[y, is] =sort (f) ;
t=tm(is(1));
