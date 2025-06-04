function B = resamph2(A,k,h,l)
%RESAMPH2     Resample and second order hold two-sided polynomial
%
% For two-sided polynomial F, the command
%  G = RESAMPH2(F,K,H,L)  returns two-sided polynomial G,
% the result of second order holding with interval L and
% resampling with period K and phase L.
%
% For more details, see FRAC/RESAMPH2, POL/RESAMPH2.

%     Author: J.Jezek, 03-Feb-2001
%     Copyright(c) 2001 by Polyx, Ltd.

ni = nargin;
if ni<1,
   error('Not enough input arguments.');
end;
if ~isa(A,'tsp'),
   error('Some argument but not 1st is invalidly two-sided polynomial.');
end;

if ni<2 | isempty(k),
   k = 2;
else
   if ~isa(k,'double') | length(k)~=1 | ~isreal(k) | ...
         floor(k)~=k | k<=0,
      error('Invalid resampling period.');
   end;
end;

if ni<3 | isempty(h),
   h = 0;
else
   if ~isa(h,'double'),
      error('Invalid resampling phase.');
   end;
end;

if ni<4 | isempty(l),
   l = k;
else
   if ~isa(l,'double'),
      error('Invalid holding interval.');
   end;
end;

oo = floor(A.o/k); o = k*oo;
Ap = shift(A.p,A.o-o,A.v);
BB = 0;
eval('BB = resamph2(Ap,k,h,l);','error(peel(lasterr));');

lh = length(h);
if lh==1,
   B = tsp(BB); B.h = A.h;
   if ~isempty(B.h), B.h = B.h*k;
   end;
   B = shift(B,oo);
else
   B = cell(size(h));
   for i = 1:lh,
      Q = tsp(BB{i}); Q.h = A.h;
      if ~isempty(Q.h), Q.h = Q.h*k;
      end;
      B{i} = shift(Q,oo);
   end;
end;

%end .. @tsp/resamph2
