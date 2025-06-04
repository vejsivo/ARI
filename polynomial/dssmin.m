function [a,b,c,d,e] = dssmin(A,B,C,D,E,tol)
%DSSMIN  Minimize the dimension of the pseudo state of a descriptor system
%
% The command
%    [a,b,c,d,e] = DSSMIN(A,B,C,D,E[,TOL])
% serves to minimize the dimension of the pseudo state of the
% descriptor system
%    Edx/dt = Ax + Bu
%       y   = Cx + Du
% Note: Any finite non-controllable or non-observable modes of the
% system are preserved.
%
% The optional parameter TOL specifies the tolerance that is used
% to separate the finite from the infinite poles. The default value
% value is TOL = 1e-8.

% Authors: H. Kwakernaak, R.C.W. Strijbos, November 27, 1998
% Copyright 1998 by PolyX Ltd.


switch nargin
case {5,6}
   if nargin == 5
      tol = 10^-8;
   else
      if ~isa(tol,'double') | length(tol)>1 
         error('Invalid tolerance.')
      end
   end
otherwise
   error('Not enough input arguments.');
end

eval('[At,Bt,Ct,Dt]=dss2ss(A,B,C,D,E,tol);', ...
   'error(peel(lasterr));');
[a,b,c,d,e]=ss2dss(At,Bt,Ct,Dt,tol);

%end .. dssmin
