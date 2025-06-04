function R = transpose(P)
%TRANSPOSE(.')    Transpose fraction
%     R = P.'  or  R = TRANSPOSE(P)
%
% For fraction P (matrix-den fraction or scalar-den fraction),
% the command returns fraction R whose (i,j)-th entry
% is  P(j,i).

%       Author:  J. Jezek  26-Jan-2000
%       Copyright (c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 25-Jul-2002 $
%                     $ Date 30-Sep-2002 $
%                     $ Date 14-Oct-2002 $

global PGLOBAL;

R = P; R.num = (P.num).'; R.den = (P.den).';
R.s = fliplr(P.s);

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

%end .. @frac/transpose
