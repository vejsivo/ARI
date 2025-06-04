function t = isreal(P,tol)
%ISREAL    Test if matrix-den fraction is real
%
% For matrix-den fraction P, the command  T = ISREAL(P)  returns 1
% if all entries of P are real, otherwise 0.
%
% If all coefficients in polynomials in P are real, then the
% result is 1, of course. However, such result can occur also
% with complex coefficients, due to a possible complex common
% factor. So, the computation is not trivial, it includes
% polynomial operations. This is why a zeroing tolerance can be
% specified by an optional input argument TOL.
%
% See also FRAC/REAL, FRAC/IMAG.

%      Author:  J. Jezek  27-Dec-1999
%      Copyright (c) 1999 by Polyx, Ltd.
%      $ Revision $     $ Date 26-Apr-2000 $
%                       $ Date 14-Oct-2002 $

global PGLOBAL;

if nargin==1,
   tol = PGLOBAL.ZEROING;
else
   if ~isa(tol,'double') | length(tol)~=1 | ...
         ~isreal(tol) | tol<0 | tol>1,
      error('Invalid tolerance.');
   end;
end;

t = logical(1);
if isreal(P.frac.num) & isreal(P.frac.den),
   return;
end;

n = P.frac.s(1)*P.frac.s(2);
for i = 1:n,
   Num = P.frac.num(i); A = real(Num); B = imag(Num);
   Den = P.frac.den(i); C = real(Den); D = imag(Den);
   if (rank([A,B;C,D],tol)>1),
      t = logical(0); return;
   end;
end;

%end .. @mdf/isreal
