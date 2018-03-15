function [A,B]=pois(nosmod,alp,beta,R,D0,D1,D2,D4);
%
% Function to create Orr-Sommerfeld matrices using Chebyshev
% pseudospectral discretization for plane Poiseuille flow
% profile
%
% nosmod = number of modes
% alp = alpha
% beta = beta
% R = Reynolds number
% D0 = zero'th derivative matrix
% D1 = first derivative matrix,
% D2 = second derivative matrix
% D3 = third derivative matrix
% D4 = fourth derivative matrix
zi=sqrt(-1);
% mean velocity
ak2=alp^2+beta^2;
Nos=nosmod+1;
Nsq=nosmod+1 ;
vec=(0:1:nosmod)';
u=(ones(length(vec),1)-cos(pi*vec/nosmod).^2);
du=-2*cos(pi*vec/nosmod);
% set up Orr-Sommerfeld matrix
B11=D2-ak2*D0;
A11=-(D4-2*ak2*D2+(ak2^2)*D0)/(zi*R);
A11=A11+alp*(u*ones(1,length(u))).*B11+alp*2*D0;
er=-200*zi;
A11=[er*D0(1,:); er*D1(1,:); A11(3:Nos-2,:);
er*D1(Nos,:); er*D0(Nos,:) ];
B11=[D0(1,:); D1(1,:); B11(3:Nos-2,:);...
D1(Nos,:); D0(Nos,:)];
% set up Squire matrix and cross-term matrix
A21=beta*(du*ones(1,length(u))).*D0(1:Nos,:);
A22=alp*(u*ones(1,length(u))).*D0-(D2-ak2*D0)/(zi*R);
B22=D0;
A22=[er*D0(1,:); A22(2:Nsq-1,:); er*D0(Nsq,:)];
A21=[zeros(1,Nos); A21(2:Nsq-1,:); zeros(1,Nos)];
% combine all the blocks
A=[A11 zeros(Nos,Nsq); A21 A22];
B=[B11 zeros(Nos,Nsq); zeros(Nsq,Nos) B22];
