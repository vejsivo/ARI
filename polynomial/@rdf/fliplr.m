function Af = fliplr(A);
%FLIPLR  Flip right-den fraction in left/right direction
%
% The command
%    FLIPLR(A) 
% returns A with its rows preserved and its columns flipped in the 
% left/right direction.              
%
% See also RDF/FLIPUD, RDF/ROT90.

%       Author: J. Jezek  07-Dec-1999
%       Copyright (c) 1999 by Polyx, Ltd.
%       $ Revision $  $ Date 25-Apr-2000 $
%                     $ Date 30-Sep-2002 $
%                     $ Date 14-Oct-2002 $

global PGLOBAL;

Af = rdf(A.frac.num,flipud(A.frac.den));

props(Af,A.frac.p,A.frac.tp);
if strcmp(PGLOBAL.COPRIME,'cop'), Af = coprime(Af);
end;
if strcmp(PGLOBAL.REDUCE,'red'), Af = reduce(Af);
else Af = smreduce(Af);
end;
if strcmp(PGLOBAL.DEFRACT,'defr'), Af = defract(Af);
end;

%end .. @rdf/fliplr
