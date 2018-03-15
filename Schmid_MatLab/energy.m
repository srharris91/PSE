function M=energy(Nos,Nsq,ak2);
%
% Program to compute the energy weight matrix for three-
% dimensional Poiseuille and Couette flows
%
% INPUT
% Nos Number of normal velocity modes
% Nsq Number of normal vorticity modes
%
% OUTPUT
% M weight matrix
M=eye(Nos+Nsq,Nos+Nsq);
Cos=two(Nos);
Dos=deven(Nos);
Wos=Dos'*Cos*Dos+ak2*Cos;
Wsq=two(Nsq);
[u,s,v]=svd(Wos); s=sqrt(diag(s));
Mos=diag(s)*u' ;
[u,s,v]=svd(Wsq); s=sqrt(diag(s));
Msq=diag (s) *u' ;
M=[Mos zeros(Nos,Nsq); zeros(Nsq,Nos) Msq];
