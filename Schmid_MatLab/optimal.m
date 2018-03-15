function [flowin,flowot,gg]=optimal(d,T,M,ak2,iflag);
%
% This function computes the initial flow structure which
% achieves the maximum transient energy growth
%
% INPUT
% d 3D Orr-Sommerfeld operator
% T time
% M energy weight matrix
% ak2 alpha^2+beta^2
% iflag flag
%       iflag = 1:
%           compute the maximum growth and
%           initial condition in time
%           interval [0, T]
%       iflag = 2:
%           compute the initial disturbance
%           yielding maximum growth at time T
%
% OUTPUT
% flowin   coefficients of optimal disturbance
%   flowin(1:Nos)           normal velocity coefficients
%   flowin(Nos+1:Nos+Nsq)   normal vorticity coefficients
%
%
%
%
% flowot coefficients of field at optimal time
%   flowot(1:Nos)           normal velocity coefficients
%   flowot(Nos+1:Nos+Nsq)   normal vorticity coefficients
%
%
%
%
global qb;
% Phase 1: Compute eigenvalues and eigenfunctions of
% Orr-Sommerfeld matrix and sort in order of descending
% imaginary part. The function nlize normalizes the
% eigenfunctions with respect to the weight matrix M.
[xs,es]=iord2(d);
xs=nlize(xs,M);
% Phase 2: Choose which modes are to be used in optimal
% calculation. Modes with imaginary part> 1 are neglected.
% Modes with imaginary part < imin are neglected as well.
ishift=1 ;
imin=-1.5;
while imag(es(ishift))>1, ishift=ishift+1; end;
[n1,n2]=n9(es,imin);
cols=(ishift:n2);
xu=xs(:,cols);
eu=es(cols);
ncols=length(cols);
fprintf('Number of modes used: %1.0f \n' ,ncols);
% Phase 3: Compute the reduced Orr-Sommerfeld operator
[qb,invF]=qbmat(M,xu,eu);
% Phase 4: Compute the time for the maximum growth using
% the built-in Matlab routine 'fmin'
if iflag==1,
    gcheck=maxer(1/100);
    gcheck=gcheck^2;
    if gcheck<1,
        tformax=O;
        mgrowth=1;
    else
         ts=T(1) ;
         tf=T(2);
         options=optimset('MaxFunEvals',500,'MaxIter',500,'TolX',1E-3');%[0 1e-3 1e-3];
         tformax=fmy(ts,tf,options);
         mgrowth=maxer(tformax);
         mgrowth=mgrowth^2;
     end;
fprintf('Time for maximum growth: %e \n',tformax);
else
    tformax=T;
end;
% Phase 5: Compute the initial condition that yields the
% maximum growth. This is obtained by
% (1) computing the matrix exponential evaluated at the
% optimal time;
% (2) computing the SVD of the matrix exponential
% exp(-i*A*t)=USV.
% The initial condition that yields the maximum growth is
% the first column of V. To convert the initial condition
% to a vector of coefficients in the eigenfunction basis
% multiply by the matrix of eigenfunctions and inv(F)
evol=expm(tformax*qb);
[U,S,V]=svd(evol);
mgrowth=S(1,1)^2;
fprintf('Maximum growth in energy: %e \n',mgrowth);
flowin=sqrt(2*ak2)*xu*invF*V(:,1);
flowot=sqrt(2*ak2)*xu*invF*U(:,1);
for i=1:100,
    tid = ts + (tf-ts)/99*(i-1);
    gg(i,2) = norm(expm(tid*qb))^2;
    gg(i,1) = tid;
end
