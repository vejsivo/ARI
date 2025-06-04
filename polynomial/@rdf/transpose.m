function G = transpose(F);
%TRANSPOSE   Transpose right-den fraction
%
% The command  G = F.'  or  G = transpose(F)  returns the transpose
% of right-denominator fraction F. The result is a left-denominator
% fraction.
%
% See also RDF/CTRANSPOSE.

%        Author:  J. Jezek  11-Nov-1999
%        Copyright(c) 1999 by Polyx, Ltd.
%        $ Revision $  $ Date 25-Apr-2000 $
%                      $ Date 30-Sep-2002 $
%                      $ Date 14-Oct-2002 $

global PGLOBAL;

G = ldf((F.frac.den).',(F.frac.num).');

props(G,F.frac.c,F.frac.tc,F.frac.r,F.frac.p,F.frac.tp);
if strcmp(PGLOBAL.COPRIME,'cop'), G = coprime(G);
end;
if strcmp(PGLOBAL.REDUCE,'red'), G = reduce(G);
else G = smreduce(G);
end;
if strcmp(PGLOBAL.DEFRACT,'defr'), G = defract(G);
end;

%end .. @rdf/transpose
