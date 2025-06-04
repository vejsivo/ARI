function [th,h,P,Q,R] = testhp3(P,Q,R,v);
%TESTHP3   Test sampling periods of three polynomials
%            [TH,H,P,Q,R] = TESTHP3(P,Q,R,V)
%
% For polynomials P,Q,R, and for V, the resulting variable 
% of some operation, the command tests whether the sampling
% periods are consistent. The resulting sampling period
% is returned in H.
%
% If all the periods are the same (up to empty ones), result
% TH = 1 . Otherwise, result TH = 0 and H is taken from V.

%      Author:  J. Jezek  24-May-2000
%      Copyright(c) 2000 by Polyx, Ltd.
%      $ Revision $  $ Date 15-Jun-2000 $

if isempty(v), vh = [];
elseif strcmp(v,'s') | strcmp(v,'p'), vh = 0;
else vh = NaN;
end;

th = 1;
if isnan(P.h) | isnan(Q.h) | isnan(R.h),
   th = 1; h = vh;
elseif isempty(P.h),
   if isempty(Q.h), h = R.h;
   elseif isempty(R.h), h = Q.h;
   elseif Q.h==R.h, h = Q.h;
   else th = 0; h = vh;
      P.h = vh; Q.h = vh; R.h = vh;
   end;
elseif isempty(Q.h),
   if isempty(R.h), h = P.h;
   elseif P.h==R.h, h = P.h;
   else th = 0; h = vh;
      P.h = vh; Q.h = vh; R.h = vh;
   end;
elseif isempty(R.h),
   if P.h==Q.h, h = P.h;
   else th = 0; h = vh;
      P.h = vh; Q.h = vh; R.h = vh;
   end;
elseif P.h==Q.h & P.h==R.h, h = P.h;
else th = 0; h = vh;
   P.h = vh; Q.h = vh; R.h = vh;
end;

%end .. private/testhp3
