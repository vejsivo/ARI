function G = mirror(F)
%MIRROR     Mirror image of fraction
%
%  G = MIRROR(F)  where F is fraction in s,p,z or z^-1,
% returns the mirror image of F:
%      G(s) = F(-s)   or    B(z) = F(z^-1)
%
% See also POL/MIRROR.

%       Author:  J. Jezek, 24-Feb-2003
%       Copyright(c) 2003 by Polyx, Ltd.

global PGLOBAL;

G = F;
eval('G.num = mirror(F.num); G.den = mirror(F.den);', ...
   'error(peel(lasterr));');
if ~isempty(G.num.v), G.v = G.num.v;
else G.v = G.den.v;
end;
G.r = 'red?';

if strcmp(PGLOBAL.COPRIME,'cop'),
   G = coprime(G);
end;
if strcmp(PGLOBAL.REDUCE,'red'),
   G = reduce(G);
else
   G = smreduce(G);
end;
if strcmp(PGLOBAL.DEFRACT,'defr'),
   G = defract(G);
end;

%end .. @frac/mirror
