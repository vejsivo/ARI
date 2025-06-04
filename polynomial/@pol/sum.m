function A = sum(P,dim,tol);
%SUM   Element-wise sum of polynomial
% 
% For a polynomial vector P (i.e., a polynomial matrix with one row
% or one column)
%    SUM(P) 
% is the sum of the elements of P. 
%
% For other polynomial matrices, SUM(P) is a row vector with the sum 
% over each column. 
%
% SUM(X,DIM) sums along the dimension DIM, where DIM is 1 or 2.
%
% SUM(X,DIM,TOL) or SUM(X,[],TOL) works with zeroing specified by 
% the input tolerance TOL.

%	Author(s): M. Hromcik, M. Sebek 31-3-98
%	Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date:  7-May-1998 12:15:34   $
%       $Revision: 3.0 $  $Date: 13-Oct-1999 12:00:00   J. Jezek  $
%                         $Date: 26-May-2000 15:00:00   J. Jezek  $
%                         $Date: 28-Feb-2003            J. Jezek  $

% Effect on other properties: 
% A.u: UserData are deleted.
% A.version: set 3.0

global PGLOBAL;

eval('P = pol(P);','error(peel(lasterr));');

ni = nargin;
if ni==1,
  dim = []; tol = PGLOBAL.ZEROING;
elseif ni==2 | isempty(tol),
  tol = PGLOBAL.ZEROING;
else
   if ~isa(tol,'double') | length(tol)~=1 | ...
         ~isreal(tol) | tol<0 | tol>1,
      error('Invalid tolerance.');
   end;
end;

Pc = P.c; 
Psi = P.s; 
Pd = P.d;

if isempty(dim),
   if isempty(Pd),
      A = pol(sum(zeros(Psi)));
      return;
   end;
   dim = find(Psi~=1);
   if isempty(dim),
      dim = 1;
   else
      dim = dim(1);
   end;
end;
if ~isa(dim,'double') | length(dim)~=1 |(dim~=1 & dim~=2),
   error('Invalid dimension.');
end;

if isempty(Pd) | isinf(Pd),
  A = pol(sum(zeros(Psi),dim));
  return;
end;

Ac = sum(P.c,dim);
A.d = Pd;
A.s = [size(Ac,1) size(Ac,2)];
A.c = Ac;

A.v = P.v;
A.h = P.h;
A.u = P.u;
A.version = 3.0;
A = class(A,'pol');

% zeroing:
me = min(abs(nonzeros(Pc)));
if ~isempty(me),
   A = pzer(A,tol*me);
end; 

%end .. @pol/sum
