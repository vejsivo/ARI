function t = isreal(F,tol)
%ISREAL    Test if right-den fraction is real
%
% For a right-denominator fraction F, the command  ISREAL(F)
% returns scalar 1 if F is real, otherwise 0.
%
% When all coefficients of numerator and of denominator of F
% are real, the result is of course 1. However, this result may
% happen even when some coefficients are complex; this effect is
% due to a possible complex common factor. So, in the case of
% complex coefficients, the computation is nontrivial; it includes
% solving a polynomial equation. It uses either the global
% zeroing tolerance or the tolerance TOL, given as an optional
% input argument.
%
% See also RDF/REAL, RDF/IMAG.

%        Author:  J. Jezek  28-Nov-1999
%        Copyright(c) 1999 by Polyx, Ltd.
%        $ Revision $  $ Date 25-Apr-2000 $
%                      $ Date 14-Oct-2002 $

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
   FF = [real(F.frac.num),imag(F.frac.num); ...
         real(F.frac.den),imag(F.frac.den)];   
   eval('t = (rank(FF,tol)==size(F,2));', ...
      'error(peel(lasterr));');
end;

%end .. @rdf/isreal
