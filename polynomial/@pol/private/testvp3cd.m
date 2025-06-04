function [tv,v,P,Q,R] = testvp3cd(P,Q,R);
%TESTPV3CD   Test variables of three polynomials,
%           continuous-time versus discrete-time
%           [TV,V,P,Q,R] = TESTPV3CD(P,Q,R)
%
% For polynomials P,Q,R, the command tests whether the variable
% symbols are consistent. The resulting variable symbol is
% returned in V.
%
% If all the symbols are the same (up to empty symbols), result
% TV = 1. Otherwise: if the symbols are of the same type
% (continuous-time or discrete-time), result TV = 0,
% else TV = -1 .

%      Author:  J. Jezek  26-May-2000
%      Copyright(c) 2000 by Polyx, Ltd.

tv = 1;
if strcmp(P.v,Q.v) & strcmp(P.v,R.v)
   v = P.v;
elseif isempty(P.v),
   if strcmp(Q.v,R.v), v = Q.v;
   elseif isempty(Q.v), v = R.v;
   elseif isempty(R.v), v = Q.v;
   else tv = 0; v = Q.v;
   end;
elseif isempty(Q.v),
   if strcmp(P.v,R.v), v = P.v;
   elseif isempty(R.v), v = P.v;
   else tv = 0; v = P.v;
   end;
elseif isempty(R.v),
   if strcmp(P.v,Q.v), v = P.v;
   else tv = 0; v = P.v;
   end;
else
   tv = 0; v = P.v;
end;

if tv==0,
   cont = {'s';'p'};
   disc = {'z';'z^-1';'q';'d'};
   I = strmatch(P.v,cont,'exact');
   J = strmatch(Q.v,cont,'exact');
   K = strmatch(R.v,cont,'exact');
   L = strmatch(P.v,disc,'exact');
   M = strmatch(Q.v,disc,'exact');
   N = strmatch(R.v,disc,'exact');
   if (~isempty(I) | ~isempty(J) | ~isempty(K)) & ...
      (~isempty(L) | ~isempty(M) | ~isempty(N)),
      tv = -1;
   end;
   P.v = v; Q.v = v; R.v = v;
end;

%end .. @pol/private/testvp3cd
