function T = tsp(R)
%TSP    Convert fraction to two-sided polynomial
%
% The command  T = TSP(R)  converts fraction R
% (matrix-den fraction or scalar-den fraction)
% to two-sided polynomial. If not possible, error.
%
% The variable symbol of R should be 'z' or 'z^-1'. If not,
% a warning is issued and 'q', 's' or 'p' is changed to 'z';
% 'd' is changed to 'z^-1'.
%
% See also FRAC/POL, FRAC/DECLASS.

%      Author:  J. Jezek  28-Jan-2000
%      Copyright(c) 2000 by Polyx, Ltd.
%      $ Revision $  $ Date 22-Jul-2002 $
%                    $ Date 14=Oct-2002 $

if ~isempty(R),
   Rc = class(R);
   if strcmp(Rc,'frac'),
      error('Argument is not convertible to tsp.');
   end;
   
   R = coprime(R);
   if strcmp(R.v,'z^-1'),
      R = reverse(R);
   elseif strcmp(R.v,'d'),
      warning('Variable symbol changed to ''z^-1''.');
      R.num.v ='z^-1'; R.den.v = 'z^-1'; R.v = 'z^-1';
      R = reverse(R);
   elseif ~isempty(R.v) & ~strcmp(R.v,'z'),
      warning('Variable symbol changed to ''z''.');
      R.num.v = 'z'; R.den.v = 'z'; R.v = 'z';
   end;

   [Dd,Dlc] = deg(R.den,'ent');
   if all(all(tdeg(R.den,'ent')==Dd)),
      P = R.num./Dlc;
      T = tsp(shift(P,-Dd,R.v));
      T.h = R.h;
   else
      error('Argument is not convertible to tsp.');
   end;
else
   T = tsp(zeros(R.s));
end;

%end .. @frac/tsp

   