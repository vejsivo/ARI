function [th,h,P,Q] = testht(P,Q);
%TESTHT   Test sampling periods of two-sided polynomials
%            [TH,H,P,Q] = TESTHT(P,Q)
%
% For two-sided polynomials P,Q, the command tests whether the
% sampling periods are consistent. The resulting sampling
% period is returned in H.
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
   th = 0; h = NaN;
end;
P.h = h; Q.h = h; P.p.h = h; Q.p.h = h;

%end .. @tsp/private/testht
