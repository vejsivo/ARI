function [tv,v,F,G] = testvf(F,G);
%TESTVF   Test variables of fractions
%           [TV,V,F,G] = TESTVF(F,G)
%
% For fractions F,G, the command tests whether
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
%      Revision  $Date: 15-Mar-2002 $
%                $Date: 14-Oct-2002 $

global PGLOBAL;

tv = 1;
if strcmp(F.v,G.v),
   v = F.v;
elseif isempty(F.v) | isempty(G.v),
   v = [F.v, G.v];
elseif (strcmp(F.v,'z') & strcmp(G.v,'z^-1')) | ...
       (strcmp(F.v,'z^-1') & strcmp(G.v,'z')),
   G = reverse(G); v = F.v;
else
   v = PGLOBAL.VARIABLE;
   F.v = v; F.num.v = v; F.den.v = v;
   G.v = v; G.num.v = v; G.den.v = v;
   if strcmp(v,'s') | strcmp(v,'p'), h = 0;
   else h = NaN;
   end;
   F.h = h; F.n.h = h; F.d.h = h;
   G.h = h; G.n.h = h; G.d.h = h;
   tv = 0;
end;

%end .. @frac/private/testvf
