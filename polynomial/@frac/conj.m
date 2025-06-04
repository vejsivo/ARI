function R = conj(P)
%CONJ    Complex conjugate of fraction
%     R = CONJ(P)
%
% For fraction P, the command returns
% fraction R whose (i,j)-th entry is
% complex conjugate of P(i,j).

%       Author:  J. Jezek  26-Jan-2000
%       Copyright (c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 30-Sep-2002 $
%                     $ Date 14-Oct-2002 $

global PGLOBAL;

R = P; R.u = [];
R.num = conj(R.num); R.den = conj(R.den);

props(R,P.p,P.tp);
Pcl = class(P);
if ~strcmp(Pcl,'frac'),
   if strcmp(PGLOBAL.COPRIME,'cop'),
      R = coprime(R);
   end;
   if strcmp(PGLOBAL.REDUCE,'red'),
      R = reduce(R);
   else
      R = smreduce(R);
   end;
   if strcmp(PGLOBAL.DEFRACT,'defr'),
      R = defract(R);
   end;
end;

%end .. @frac/conj
