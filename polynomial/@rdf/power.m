function An = power(A,n,tol);
%POWER (.^) Element-wise power of right-den fraction
%           A.^n
%
% AN = A.^N or AN = POWER(A,N) denotes element-by-element powers.
% The (i,j)-the element of the result is A(i,j)^N(i,j). The right
% fraction A and the integer matrix N must have the same sizes
% unless one of them is scalar. The scalar operates into every
% element of the other matrix.
%
% If any element of matrix N is negative then the corresponding
% element of A must be nonzero. Note that the individual elements
% A(i,j) are not apparent in the form A=NUM/DEN.
%
% AN = POWER(A,N,TOL) works with zeroing specified by the input
% relative tolerance TOL.
%
% See also FRAC/MPOWER, RDF/TIMES.

%       Author:  J. Jezek  28-Jan-2000
%       Copyright(c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 25-Apr-2000 $

global PGLOBAL;

na  = nargin;
if na<2,
   error('Not enough input arguments.');
elseif na==2 | isempty(tol),
   tol = PGLOBAL.ZEROING;
else
   if ~isa(tol,'double'),
      error('Invalid tolerance.');
   end;
end;   

if ~isa(n,'double'),
   error('Invalid power; must be integer.');
end;

cop = PGLOBAL.COPRIME; PGLOBAL.COPRIME = 'cop';
eval('An = rdf(power(mdf(A),n,tol));','error(peel(lasterr));');
PGLOBAL.COPRIME = cop;

%end .. @rdf/power
