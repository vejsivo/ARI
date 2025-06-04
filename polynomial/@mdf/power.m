function An = power(A,n,tol);
%POWER  (.^) Element-wise power of matrix-den fraction
%           A.^n
%
% AN = A.^N or AN = POWER(A,N) denotes element-by-element powers. 
% The matrix-den fraction A and the integer matrix N must have the
% same sizes unless one of them is a scalar. The scalar operates
% with every element of the other matrix.
%
% If any element of matrix N is negative then the corresponding
% element of A must be nonzero.
%
% AN = POWER(A,N,TOL) works with zeroing specified by the input
% relative tolerance TOL.
%
% See also FRAC/MPOWER, MDF/TIMES.

%       Author:  J. Jezek  06-Jan-2000
%       Copyright(c) by Polyx, Ltd.
%       $ Revision $  $ Date 26-Apr-2000 $
%                     $ Date 30-Sep-2002 $
%                     $ Date 14-Oct-2002 $

global PGLOBAL;

na  = nargin;
if na<2,
   error('Not enough input arguments.');
elseif na==2 | isempty(tol),
   tol = PGLOBAL.ZEROING;
else
   if ~isa(tol,'double') | length(tol)~=1 | ...
         ~isreal(tol) | tol<0 | tol>1,
      error('Invalid tolerance.');
   end;
end;   

if ~isa(n,'double') | ndims(n)>2 | ~isreal(n) |  ...
      (~isempty(n) & any(any(floor(n)~=n)) ),
   error('Invalid power; must be integer.');
end;

As = A.frac.s; ns = size(n);
if any(As~=ns),
   if all(As==1),
      A = repmat(A,ns);
   elseif all(ns==1),
      n = repmat(n,As);
   else
      error('Matrices not of the same dimensions.');
   end;
end;

allnpos = logical(1);
As1 = A.frac.s(1); As2 = A.frac.s(2);
for i = 1:As1,
   for j = 1:As2,
      if n(i,j)<0,
         allnpos = logical(0);
         if eq(A.frac.num(i,j),0,tol),
            error('Zero powered to a negative power.');
         end;
         B = A.frac.num(i,j); A.frac.num(i,j) = A.frac.den(i,j);
         A.frac.den(i,j) = B;
         n(i,j) = abs(n(i,j));
      end;
   end;
end;

An = mdf(power(A.frac.num,n,tol),power(A.frac.den,n,tol));

isproper(An);
if strcmp(A.frac.c,'cop'),
   props(An,'cop',A.frac.tc);
end;
if strcmp(A.frac.r,'red') & allnpos,
   props(An,'red');
end;

if strcmp(PGLOBAL.COPRIME,'cop'),
   An = coprime(An,tol);
end;
if strcmp(PGLOBAL.REDUCE,'red'),
   An = reduce(An);
else
   An = smreduce(An);
end;
if strcmp(PGLOBAL.DEFRACT,'defr'),
   An = defract(An);
end;

%end .. @mdf/power
