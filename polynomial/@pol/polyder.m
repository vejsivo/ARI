function P = polyder(P, N)
%POLYDER  Derivative of a polynomial matrix
%
% POLYDER(P) computes the first derivative of the polynomial matrix P.
% POLYDER(P,N) computes the Nth derivative of the polynomial matrix P.

%    Author: D. Henrion, September 20, 1998.
%    Modified by D. Henrion, September 10, 1999.
%    Updated to 3.0 by D. Henrion, August 11, 2000.
%    Copyright 1998-2000 by Polyx, Ltd.
%    Modified by J. Jezek, September 04, 2001.

if nargin < 2,
 N = 1;
elseif ~isa(N, 'double'),
 error('Invalid 2nd argument.');
elseif length(N) > 1 | N < 0,
 error('Invalid 2nd argument; must be a positive scalar.');
end;

P = pol(P);
degP = deg(P);
[rP, cP] = size(P);

if N > 0,

 if isempty(degP) | degP <= 0,

  P = pol(zeros(rP, cP));

 else

  Pc = zeros(rP, cP*degP);

  for i = 1:degP,
   Pc(:, 1+(i-1)*cP:i*cP) = i * P.c(:,:,i+1);
  end;

  P = pol(Pc, degP - 1, P.v);

  P = polyder(P, N - 1);

 end;

end;

%end .. @pol/polyder
