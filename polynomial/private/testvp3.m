function [tv,v,P,Q,R] = testvp3(P,Q,R);
%TESTVP3   Test variables of three polynomials
%           [TV,V,P,Q,R] = TESTVP3(P,Q,R)
%
% For polynomials P,Q,R, the command tests whether the variable
% symbols are consistent. The resulting variable symbol is
% returned in V.
%
% If all the symbols are the same (up to empty symbols), result
% TV = 1. If the symbols are mixture of z, z^-1 and possible
% empty, result TV = 2 and  V = PGLOBAL.DISCRVAR .
% Otherwise, result TV = 0 and  V = PGLOBAL.VARIABLE .

%      Author:  J. Jezek  24-May-2000
%      Copyright(c) 2000 by Polyx, Ltd.

global PGLOBAL;

tv = 1;
if strcmp(P.v,Q.v) & strcmp(P.v,R.v)
   v = P.v;
elseif isempty(P.v),
   if strcmp(Q.v,R.v), v = Q.v;
   elseif isempty(Q.v), v = R.v;
   elseif isempty(R.v), v = Q.v;
   else tv = 0;
   end;
elseif isempty(Q.v),
   if strcmp(P.v,R.v), v = P.v;
   elseif isempty(R.v), v = P.v;
   else tv = 0;
   end;
elseif isempty(R.v),
   if strcmp(P.v,Q.v), v = P.v;
   else tv = 0;
   end;
else
   tv = 0;
end;

if tv==0,
   if (isempty(P.v) | strcmp(P.v,'z') | strcmp(P.v,'z^-1')) & ...
      (isempty(Q.v) | strcmp(Q.v,'z') | strcmp(Q.v,'z^-1')) & ...
      (isempty(R.v) | strcmp(R.v,'z') | strcmp(R.v,'z^-1')),
         tv = 2; v = PGLOBAL.DISCRVAR;
   else      
       v = PGLOBAL.VARIABLE;
      if strcmp(v,'s') | strcmp(v,'p'), h = 0;
      else h = NaN;
      end;
      if isempty(P.d) | P.d<=0, P.v = ''; P.h = [];
      else P.v = v; P.h = h;
      end;
      if isempty(Q.d) | Q.d<=0, Q.v = ''; Q.h = [];
      else Q.v = v; Q.h = h;
      end;
      if isempty(R.d) | R.d<=0, R.v = ''; R.h = [];
      else R.v = v; R.h = h;
      end;
   end;
end;      

%end .. private/testvp3
