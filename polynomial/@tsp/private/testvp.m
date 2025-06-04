function [tv,v,P,Q] = testvp(P,Q);
%TESTVP   Test variables of polynomials
%           [TV,V,P,Q] = TESTVP(P,Q)
%
% For polynomials P,Q, the command tests whether the variable
% symbols are consistent. The resulting variable symbol is
% returned in V.
%
% If the symbols are the same or some of them is empty, result TV = 1.
% If one of the symbols is 'z' and the other 'z^-1', result TV = 2
% and  V = PGLOBAL.DISCRVAR . Otherwise, result TV = 0 and
% V = PGLOBAL.VARIABLE .

%      Author:  J. Jezek  12-May-2000
%      Copyright(c) 2000 by Polyx, Ltd.

global PGLOBAL;

if strcmp(P.v,Q.v),
   tv = 1; v = P.v;
elseif isempty(P.v) | isempty (Q.v),
   tv = 1; v = [P.v,Q.v];
elseif (strcmp(P.v,'z') & strcmp(Q.v,'z^-1')) | ...
       (strcmp(P.v,'z^-1') & strcmp(Q.v,'z')),
   tv = 2; v = PGLOBAL.DISCRVAR;
else
   tv = 0; v = PGLOBAL.VARIABLE;
   if strcmp(v,'s') | strcmp(v,'p'), h = 0;
   else h = NaN;
   end;
   if isempty(P.d) | P.d<=0, P.v = ''; P.h = [];
   else P.v = v; P.h = h;
   end;
   if isempty(Q.d) | Q.d<=0, Q.v = ''; Q.h = [];
   else Q.v = v; Q.h = h;
   end;
end;

%end .. @tsp/private/testvp
