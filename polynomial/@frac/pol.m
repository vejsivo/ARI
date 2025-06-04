function P = pol(F)
%POL    Convert fraction to polynomial
%
% The command  P = POL(F)  converts fraction F to polynomial.
% If not possible, error.
%
% See also FRAC/DECLASS.

%        Author:  J. Jezek  10-Jan-2000
%        Copyright(c) 2000 by Polyx, Ltd.
%        $ Revision $  $ Date 22-Jul-2002 $
%                      $ Date 14-Oct-2002 $

if ~isempty(F),
   Fc = class(F);
   if strcmp(Fc,'frac'),
      error('Argument is not convertible to polynomial.');
   end;
   F = reduce(coprime(F));
   Dd = deg(F.den);
   if Dd==0,
      P = F.num;
   elseif (strcmp(F.v,'z') | strcmp(F.v,'z^-1')) & ...
         tdeg(F.den)==Dd & Dd>=deg(F.num),
      F = reverse(F); P = pol(F);
      return;
   else
      error('Argument is not convertible to polynomial.');
   end;
else
   P = pol(zeros(F.s));
end;

%end .. @frac/pol
