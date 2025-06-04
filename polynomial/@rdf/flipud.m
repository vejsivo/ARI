function Af = flipud(A);
%FLIPUD  Flip right-den fraction in up/down direction
%
% The command
%    FLIPUD(A) 
% returns A with its columns preserved and its rows flipped
% in the up/down direction.
%
% See also RDF/FLIPLR, RDF/ROT90.

%       Author: J. Jezek  07-Nov-1999
%       Copyright (c) 1999 by Polyx, Ltd.
%       $ Revision $  $ Date 25-Apr-2000 $
%                     $ Date 30-Sep-2002 $
%                     $ Date 14-Oct-2002 $

global PGLOBAL;

Af = rdf(flipud(A.frac.num),A.frac.den);

props(Af,A.frac.p,A.frac.tp);
if strcmp(PGLOBAL.COPRIME,'cop'), Af = coprime(Af);
end;
if strcmp(PGLOBAL.REDUCE,'red'), Af = reduce(Af);
else Af = smreduce(Af);
end;
if strcmp(PGLOBAL.DEFRACT,'defr'), Af = defract(Af);
end;

%end .. @rdf/flipud
