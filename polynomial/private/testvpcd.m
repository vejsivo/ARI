function [tv,v,P,Q] = testvpcd(P,Q);
%TESTVPCD   Test variables of polynomials,
%          continuous-time versus discrete-time
%           [TV,V,P,Q] = TESTVPCD(P,Q)
%
% For polynomials P,Q, the command tests whether the variable
% symbols are consistent. The resulting variable symbol is
% returned in V.
%
% If the symbols are the same or some of them is empty, result TV = 1.
% Otherwise: if the symbols are of the same type (continuous-time
% or discrete-time), result TV = 0, else TV = -1 .

%      Author:  J. Jezek  25-May-2000
%      Copyright(c) 2000 by Polyx, Ltd.

if strcmp(P.v,Q.v),
   tv = 1; v = P.v;
elseif isempty(P.v) | isempty (Q.v),
   tv = 1; v = [P.v,Q.v];
else
   I = strmatch(P.v,{'s';'p'},'exact');
   J = strmatch(Q.v,{'s';'p'},'exact');
   K = strmatch(P.v,{'z';'z^-1';'d';'q'},'exact');
   L = strmatch(Q.v,{'z';'z^-1';'d';'q'},'exact');
   if (~isempty(I) & ~isempty(L)) | ...
         (~isempty(J) & ~isempty(K)),
      tv = -1;
   else
      tv = 0;
   end;
   v = P.v; Q.v = v;
end;

%end .. private/testvpcd
