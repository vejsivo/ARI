function [v,h,N,D] = testvhnd(N,D)
%TESTVHND    Test variable and sampling period of polynomials
%            to be a numerator and a denominator of a fraction
%
% The command   [V,H,N,D] = TESTVHND(N,D)  tests variables
% and sampling periods of polynomials N,D, which are to be
% a numerator and a denominator of a fraction. Admissible
% is also a combination of 'z' and 'z^-1'; this case is
% converted to 'z'.

%      Author: J. Jezek, 08-Aug-2001
%      Copyright(c) 2001 by Polyx, Ltd.

[tv,v,N,D] = testvp(N,D);
if tv==2,
   [th,h,N,D] = testhp(N,D,v);
   if ~th,
      warning('Inconsistent sampling periods.');
   end;
   if strcmp(D.v,'z^-1'),
      dg = deg(D);
      N = shift(N,dg); D = rev(D,dg);
      D.v = 'z';
   else
      dg = deg(N);
      N = rev(N,dg); D = shift(D,dg);
      N.v = 'z';
   end;
   v = 'z';
elseif ~tv,
   warning('Inconsistent variables.');
end;

[th,h,N,D] = testhp(N,D,v);
if ~th,
   warning('Inconsistent sampling periods.');
end;

%end .. private/testvhnd
