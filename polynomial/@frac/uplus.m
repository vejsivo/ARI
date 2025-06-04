function R = uplus(P);
%UPLUS(+)    Unary plus of fraction
%     R = +P  or  R = UPLUS(P)
%
% For scalar fraction P, the command returns
% fraction R whose (i,j)-th entry is +P(i,j).

%     Author:  J. Jezek  26-Jan-2000
%     Copyright(c) 2000 by Polyx, Ltd.
%     $ Revision $  $ Date 26-Jul-2002 $
%                   $ Date 30-Sep-2002 $

global PGLOBAL;

R = P;
R.u = [];

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
   if strcmp(PGLOBAL.DEFRACT,'defr');
      R = defract(R);
   end;
end;

%end .. @frac/uplus
