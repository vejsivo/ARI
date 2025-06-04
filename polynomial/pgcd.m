function [G,Q] = pgcd(P,dim,tol)
%PGCD    Polynomial greatest common divisor
%
% For a polynomial vector P (i.e., a polynomial matrix with one row
% or one column), the command  G = PGCD(P)  returns the greatest common
% divisor of the elements of P. The command  [G,Q] = PLCM(P) returns
% also polynomial vector Q such that i-th entry, it holds G*Q(i)=P(i).
%
% For other polynomial matrix P, the returned G is a row vector with
% greatest common divisors over each column. The returned Q is 
% a polynomial matrix of the same size as P. For (i,j)-th entry,
% it holds  G(j)*Q(i,j)=P(i,j).
%
% PGCD(P,DIM) computes gcd along the dimension DIM,
% where DIM is 1 or 2.
%
% PGCD(P,DIM,TOL)  or  PGCD(P,[],TOL)  works with zeroing tolerance
% TOL.

% Note: This macro must NOT reside in directory @POL !

%       Author:  J. Jezek  25-Jan-2000
%       Copyright(c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 26-May-2000 $

global PGLOBAL;
eval('PGLOBAL.ZEROING;', 'painit;');

ni = nargin;
if ni==0
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
if Ps1==1 & (isempty(dim) | dim==2),
   n = Ps2;
elseif Ps2==1 & (isempty(dim) | dim==1),
   n = Ps1;
elseif isempty(dim) | dim==1,
   G = pol(ones(1,Ps2)); Q = P;
   if Ps2>0,
      for j = 1:Ps2,
         [G(j),Q(:,j)] = pgcd(P(:,j),[],tol);
      end;
   end;
   return;
else
   G = pol(ones(Ps1,1)); Q = P;
   if Ps1>0,
      for i = 1:Ps1,
         [G(i),Q(i,:)] = pgcd(P(i,:),[],tol);
      end;
   end;
   return;
end;

if n==0,
   G = pol(0);
else
   PP = cell(1,n);
   for i = 1:n,
      PP{i} = P(i);
   end;
   G = grd(PP,tol);
   G.h = P.h;
end;

if nargout==2,
   Q = pol(zeros(P.s));
   if n>0,
      for i = 1:n,
         Q(i) = rdiv(P(i),G,tol/10);
      end;
   end;
   Q.h = P.h;
end;

%end .. pgcd
