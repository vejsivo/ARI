function [th,h,F,G] = testhf(F,G,v);
%TESTHF   Test sampling periods of left-den fractions
%            [TH,H,F,G] = TESTHF(F,G,V)
%
% For left-den fractions F,G, and for V, the resulting variable 
% of the operation, the command tests whether the sampling
% periods are consistent. The resulting sampling period
% is returned in H.
%
% If the periods are the same or one of them is empty, result
% TH = 1 . Otherwise, result TH = 0 and H is taken from V.

%      Author:  J. Jezek  22-May-2000
%      Copyright(c) 2000 by Polyx, Ltd.
%      $ Revision $  $ Date 15-Jun-2000 $
%                    $ Date 14-Oct-2002 $

if isnan(F.frac.h) | isnan(G.frac.h),
   th = 1; h = NaN;
elseif isempty(F.frac.h),
   th = 1; h = G.frac.h;
elseif isempty(G.frac.h),
   th = 1; h = F.frac.h;
elseif F.frac.h == G.frac.h,
   th = 1; h = F.frac.h;
else
   th = 0;
   if isempty(v), h = [];
   elseif strcmp(v,'s') | strcmp(v,'p'), h = 0;
   else h = NaN;
   end;
end;
F.frac.h = h; F.frac.num.h = h; F.frac.den.h = h;
G.frac.h = h; G.frac.num.h = h; G.frac.den.h = h;

%end .. @ldf/private/testhf
