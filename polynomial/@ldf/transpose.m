function G = transpose(F);
%TRANSPOSE   Transpose of left-den fraction
%
% The command  G = F.'  or  G = TRANSPOSE(F)  returns the transpose
% of left-denominator fraction F. The result is a right-denominator
% fraction.
%
% See also LDF/CTRANSPOSE.

%        Author:  J. Jezek  11-Nov-1999
%        Copyright(c) 1999 by Polyx, Ltd.
%        $ Revision $  $ Date 21-Apr-2000 $
%                      $ Date 30-Sep-2002 $
%                      $ Date 14-Oct-2002 $

global PGLOBAL;

G = rdf((F.frac.num).',(F.frac.den).');

props(G,F.frac.c,F.frac.tc,F.frac.r,F.frac.p,F.frac.tp);
if strcmp(PGLOBAL.COPRIME,'cop'), G = coprime(G);
end;
if strcmp(PGLOBAL.REDUCE,'red'), G = reduce(G);
else G = smreduce(G);
end;
if strcmp(PGLOBAL.DEFRACT,'defr'), G = defract(G);
end;

%end .. @ldf/transpose
