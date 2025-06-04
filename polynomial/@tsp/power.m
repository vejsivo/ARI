function [An,C] = power(A,n,tol)
%POWER    Element-wise power of two-sided polynomial
%                   AN = A.^N
%
% AN = A.^N or AN = POWER(A,N) denotes element-by-element powers. 
% The tsp matrix A and the integer matrix N must have the same
% sizes unless one of them is a scalar. The scalar operates with 
% every element of the other matrix. 
%
% If all elements of N are nonnegative then the result is a two-
% sided polynomial. If any of them is negative then the result is
% a matrix-denominator fraction. In such a case, the corresponding
% element of A must be nonzero. The class of the result can be
% obtained in an optional output argument.
%
% AN = POWER(A,N,TOL) works with zeroing specified by the input
% relative tolerance TOL.
%
% See also TSP/MPOWER, TSP/TIMES.

%       Author:  J. Jezek  11-8-99
%       Copyright (c) 1999 by Polyx, Ltd.
%       $Revision: 3.0 $  $Date: 21-Apr-2000 $
%                         $Date: 24-May-2000 $
%                         $Date: 31-Oct-2000 $

global PGLOBAL;

na = nargin;
if na==2,
   tol = PGLOBAL.ZEROING;
elseif na==3,
   if ~isa(tol,'double'),
      error('Invalid 3rd argument.');
   end;
else
   error ('Not enough input arguments.');   
end;

if ~isa(n,'double') | ndims(n)>2 | ~isreal(n) |  ...
      (~isempty(n) & any(any(floor(n)~=n)) ),
   error('Invalid power; must be integer.');
end;

if ~all(size(A)==size(n)),
   if all(size(A)==1),
      A = repmat(A,size(n));
   elseif all(size(n)==1),
      n = repmat(n,size(A));
   else
      error('Matrices not of the same dimensions.');
   end;
end;

PP = 0;
eval('PP  = power(A.p,n,tol);','error(peel(lasterr));');
Aon = A.o*n;
if isa(PP,'pol'),
   An = tsp(PP); An.h = A.h;
   An = shift(An,Aon);
else      %isa(PP,'mdf')
   [ns1,ns2] = size(n);
   for i = 1:ns1,
      for j = 1:ns2,
         Aonij = Aon(i,j);
         if Aonij>=0,
            PP.n(i,j) = shift(PP.n(i,j),Aonij,A.v);
         else
            PP.d(i,j) = shift(PP.d(i,j),abs(Aonij),A.v);
         end;
      end;
   end;
   An = PP; An.h = A.h;
   if strcmp(PGLOBAL.COPRIME,'cop'), An = coprime(An);
   end;
   if strcmp(PGLOBAL.REDUCE,'red'), An = reduce(An);
   else An = smreduce(An);
   end;
end;
C = class(An);

%end .. @tsp/power

