function NP = npos(P)
%NPOS  Nonpositive part of a polynomial
%
% The command
%      NP = NPOS(P)
% for the polynomial P, returns the part with
% nonpositive terms only. If P is polynomial in
% other variable than z^-1, the result is only
% the absolute term, if any.
%
% See also POL/NNEG, POL/NEG, POL/POS.

%       Author:  J. Jezek  11-8-99
%       Copyright 1999 by Polyx, Ltd.

NP = P;
if ~isempty(P.v) & ~strcmp(P.v,'z^-1'),
   Ps = P.s; U = P.c;
   V = U(Ps(1),Ps(2),1); NP = pol(V);
end;

%end .. @pol/npos
