function [tv,v,F,G] = testvf(F,G);
%TESTVF   Test variables of right-den fractions
%           [TV,V,F,G] = TESTVF(F,G)
%
% For right-den fractions F,G, the command tests whether
% the variable symbols are consistent. The resulting
% variable symbol is returned in V.
%
% If the symbols are the same or one of them is empty,
% result TV = 1. When not, the symbols are changed to
% the standard one, resulting TV = 0. However, if one
% symbol is 'z' and the other 'z^-1' then they are taken
% as consistent. The resulting symbol is taken from P,
% resulting TV = 1.

%      Author:  J. Jezek  04-Feb-2000
%      Copyright(c) 2000 by Polyx, Ltd.
%      $ Revision $  $ Date 14-Oct-2002 $

global PGLOBAL;

tv = 1;
if strcmp(F.frac.v,G.frac.v),
   v = F.frac.v;
elseif isempty(F.frac.v) | isempty(G.frac.v),
   v = [F.frac.v, G.frac.v];
elseif (strcmp(F.frac.v,'z') & strcmp(G.frac.v,'z^-1')) | ...
       (strcmp(F.frac.v,'z^-1') & strcmp(G.frac.v,'z')),
   G = reverse(G); v = F.frac.v;
else
   v = PGLOBAL.VARIABLE;
   F.frac.v = v; F.frac.num.v = v; F.frac.den.v = v;
   G.frac.v = v; G.frac.num.v = v; G.frac.den.v = v;
   if strcmp(v,'s') | strcmp(v,'p'), h = 0;
   else h = NaN;
   end;
   F.frac.h = h; F.frac.num.h = h; F.frac.den.h = h;
   G.frac.h = h; G.frac.num.h = h; G.frac.den.h = h;
   tv = 0;
end;

%end .. @rdf/private/testvf
