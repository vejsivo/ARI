function [th,h,F,G] = testhf(F,G,v);
%TESTHF   Test sampling periods of fractions
%            [TH,H,F,G] = TESTHF(F,G,V)
%
% For fractions F,G, and for V, the resulting variable 
% of the operation, the command tests whether the sampling
% periods are consistent. The resulting sampling period
% is returned in H.
%
% If the periods are the same or one of them is empty, result
% TH = 1 . Otherwise, result TH = 0 and H is taken from V.

%      Author:  J. Jezek  22-May-2000
%      Copyright(c) 2000 by Polyx, Ltd.
%      $ Revision $  $ Date 15-Jun-2000 $
%                    $ Date 19-Sep-2001 $
%                    $ Date 14-Oct-2002 $

if isnan(F.h) | isnan(G.h),
   th = 1; h = NaN;
elseif isempty(F.h),
   th = 1; h = G.h;
elseif isempty(G.h),
   th = 1; h = F.h;
elseif F.h == G.h,
   th = 1; h = F.h;
else
   th = 0;
   if isempty(v), h = [];
   elseif strcmp(v,'s') | strcmp(v,'p'), h = 0;
   else h = NaN;
   end;
end;
F.h = h; F.num.h = h; F.den.h = h;
G.h = h; G.num.h = h; G.den.h = h;

%end .. @frac/private/testhf
