function T = tsp(F)
%TSP     Convert right-den fraction to two-sided polynomial
%
% The command  T = TSP(F)  converts right-den fraction F to
% two-sided polynomial. If not possible, error.
%
% The variable symbol of F should be 'z' or 'z^-1'. If not,
% a warning is issued and 'q','s','p' is changed to 'z',
% 'd' is changed to 'z^-1'.
%
% See also FRAC/POL, RDF/DECLASS.

%       Author:  J. Jezek  10-Jan-2000
%       Copyright(c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 28-Feb-2003 $

if ~isempty(F),
   F = reduce(coprime(F));
   if strcmp(F.frac.v,'z^-1'),
      F = reverse(F);
   elseif strcmp(F.frac.v,'d'),
      warning('Variable symbol changed to ''z^-1''.');
      F.frac.num.v ='z^-1'; F.frac.den.v = 'z^-1'; F.frac.v = 'z^-1';
      F = reverse(F);
   elseif ~isempty(F.frac.v) & ~strcmp(F.frac.v,'z'),
      warning('Variable symbol changed to ''z''.');
      F.frac.num.v = 'z'; F.frac.den.v = 'z'; F.frac.v = 'z';
   end;
   Dd = deg(F.frac.den,'row');
   if all(tdeg(F.frac.den,'row')==Dd),
      Dd = repmat(Dd,1,F.frac.s(2));
      T = tsp(shift(F.frac.num,-Dd,F.frac.v));
   else
      error('Argument is not convertible to tsp.');
   end;
else
   T = tsp(zeros(F.frac.s));
end;

%end .. @rdf/tsp
