function [M,Q] = plcm(P,dim,tol)
%PLCM    Polynomial least common multiple
%
% For a polynomial vector P (i.e., a polynomial matrix with one row
% or one column), the command  M = PLCM(P)  returns the least common
% multiple of the elements of P. The command  [M,Q] = PLCM(P) returns
% also polynomial vector Q such that i-th entry, it holds P(i)*Q(i)=M.
%
% For other polynomial matrix P, the returned M is a row vector with
% least common multiples over each column. The returned Q is 
% a polynomial matrix of the same size as P. For (i,j)-th entry,
% it holds  P(i,j)*Q(i,j)=M(j).
%
% PLCM(P,DIM) computes lcm along the dimension DIM,
% where DIM is 1 or 2.
%
% PLCM(P,DIM,TOL)  or  PLCM(P,[],TOL)  works with zeroing tolerance
% TOL.

% Note: This macro must NOT reside in directory @POL !

%       Author:  J. Jezek  05-Jan-2000
%       Copyright(c) 2000 by Polyx, Ltd.

global PGLOBAL;
eval('PGLOBAL.ZEROING;', 'painit;');

ni = nargin;
if ni==0,
   error('Not enough input arguments.');
elseif ni==1,
   dim = []; tol = PGLOBAL.ZEROING;
elseif ni==2 | isempty(tol),
   tol = PGLOBAL.ZEROING;
elseif ~isa(tol,'double') | length(tol)~=1 | ~isreal(tol) | ...
      tol<0 | tol>1,
   error('Invalid tolerance.');
end;
if ~isempty(dim) & (~isa(dim,'double') | length(dim)~=1 | ...
      (dim~=1 & dim~=2)),
   error('Invalid dimension.');
end;

eval('P = pol(P);','error(peel(lasterr));');

Ps1 = P.s(1); Ps2 = P.s(2);
%if Ps1==0 & Ps2==0,
%  M = pol(1); Q = P;
%   return;
if Ps1==1 & (isempty(dim) | dim==2),
   n = Ps2;
elseif Ps2==1 & (isempty(dim) | dim==1),
   n = Ps1;
elseif isempty(dim) | dim==1,
   M = pol(zeros(1,Ps2)); Q = P;
   if Ps2>0,
      for j = 1:Ps2,
         [M(j),Q(:,j)] = plcm(P(:,j),[],tol);
      end;
   end;
   M.h = P.h;
   return;
else
   M = pol(zeros(Ps1,1)); Q = P;
   if Ps1>0,
      for i = 1:Ps1,
         [M(i),Q(i,:)] = plcm(P(i,:),[],tol);
      end;
   end;
   M.h = P.h;
   return;
end;

if n==0,
   M = pol(1); Q = P;
   return;
end;
MM = horzcat(diag(P(1:n-1)),zeros(n-1,1)) - ...
     horzcat(zeros(n-1,1),diag(P(2:n)));
Q = null(MM,tol);
Q = reshape(Q,Ps1,Ps2);
M = times(P(1),Q(1),tol);

Md = M.d;
if ~isempty(Md) & Md>=0,
   lead = 1/M{Md};
   M = M*lead; Q = Q*lead;
end;

%end .. plcm



