function R = uminus(P);
%UMINUS(-)    Unary minus of fraction
%     R = -P  or  R = UMINUS(P)
%
% For fraction P, the command returns
% fraction R whose (i,j)-th entry is -P(i,j).

%     Author:  J. Jezek  26-Jan-2000
%     Copyright(c) 2000 by Polyx, Ltd.
%     $ Revision $  $ Date 25-Jul-2002 $
%                   $ Date 30-Sep-2002 $
%                   $ Date 14-Oct-2002 $

global PGLOBAL;

R = P; R.u = []; R.num = -R.num;

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

%end .. @frac/uminus
