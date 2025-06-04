function [th,h,P,Q] = testhp(P,Q,v);
%TESTHP   Test sampling periods of polynomials
%            [TH,H,P,Q] = TESTHP(P,Q,V)
%
% For polynomials P,Q, and for V, the resulting variable 
% of the operation, the command tests whether the sampling
% periods are consistent. The resulting sampling period
% is returned in H.
%
% If the periods are the same or one of them is empty, result
% TH = 1 . Otherwise, result TH = 0 and H is taken from V.

%      Author:  J. Jezek  22-May-2000
%      Copyright(c) 2000 by Polyx, Ltd.
%      $ Revision $  $ Date 15-Jun-2000 $ 

if isnan(P.h) | isnan(Q.h),
   th = 1; h = NaN;
elseif isempty(P.h),
   th = 1; h = Q.h;
elseif isempty(Q.h),
   th = 1; h = P.h;
elseif P.h==Q.h,
   th = 1; h = P.h;
else
   th = 0;
   if isempty(v), h = [];
   elseif strcmp(v,'s') | strcmp(v,'p'), h = 0;
   else h = NaN;
   end;
end;
if ~isempty(P.v), P.h = h;
else P.h = [];
end;
if ~isempty(Q.v), Q.h = h;   
else Q.h = [];
end;

%end .. @rdf/private/testhp
