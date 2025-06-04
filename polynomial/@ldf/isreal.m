function t = isreal(F,tol)
%ISREAL    Test if left-den fraction is real
%
% For a left-denominator fraction F, the command  ISREAL(F)
% returns scalar 1 if F is real, otherwise 0.
%
% When all coefficients of the numerator and of the denominator of
% F are real, the result is of course 1. However, this result may
% happen even when some coefficients are complex; this effect is
% due to a possible complex common factor. So, in the case of
% complex coefficients, the computation is nontrivial; it includes
% solving a polynomial equation. It uses either the global
% zeroing tolerance or the tolerance TOL, given as an optional
% input argument.
%
% See also LDF/REAL, LDF/IMAG.

%        Author:  J. Jezek  28-Nov-1999
%        Copyright(c) 1999 by Polyx, Ltd.
%        $ Revision $  $ Date 14-Oct-2002 $

global PGLOBAL;
if nargin==2,
   if ~isa(tol,'double'),
      error('Invalid tolerance.');
   end;
else
   tol = PGLOBAL.ZEROING;
end;

if isreal(F.frac.num) & isreal (F.frac.den),
   t = logical(1);
else
   FF = [real(F.frac.den),real(F.frac.num); ...
         imag(F.frac.den),imag(F.frac.num)];   
   eval('t = (rank(FF,tol)==size(F,1));', ...
      'error(peel(lasterr));');
end;

%end .. @ldf/isreal
