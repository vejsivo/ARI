function r = elliana(p0,pq)
% The function
%
%    R = ELLIANA(P0,[P1 P2 .. PN])
%
% computes the stability radius R of the continuous-time polynomial
% P = P0 + Q1*P1 + Q2*P2 + .. + QN*PN with real ellipsoidal uncertainty
% Q = [Q1 Q2 .. QN], i.e. R is the maximum Euclidean norm of vector Q
% such that P is robustly stable.
  
% Written by D. Henrion, October 15, 2002
% Last modified by D. Henrion, October 23, 2002
% Copyright 2002 by PolyX, Ltd
  
global PGLOBAL;

tol = 1e4*PGLOBAL.ZEROING;
verbose = strcmp(PGLOBAL.VERBOSE,'yes');

p0 = pol(p0);

if ~isstable(p0),
 error('Nominal polynomial is not stable');
end;

if nargin < 2,

  % no uncertainty
  r = inf;

else
  
  pq = pol(pq);
  
  if (~strcmp(symbol(p0),'s') & ~strcmp(symbol(p0),'q') & ...
    ~isempty(symbol(p0))) | (~strcmp(symbol(pq),'s') & ...
    ~strcmp(symbol(pq),'q') & ~isempty(symbol(pq))),
    error('Continuous-time polynomials only');
  end;

  spq = size(pq);
  if min(spq) > 1,
    error('Second input argument must be a polynomial vector');
  end;
  n = max(spq);

  if deg(pq) > deg(p0)
    error('Invalid degrees of input polynomials');
  end;
  
  % make nominal polynomial monic
  pq = pq / lcoef(p0);
  p0 = p0 / lcoef(p0);
  
  % extract real (x) and imaginary (y) parts
  x = pol(zeros(1,n+1)); y = pol(zeros(1,n+1));
  for i = 0:n,
    if i == 0, p = p0; else p = pol(pq(i)); end;
    dp = deg(p);
    coef = p{0:2:dp}; d = length(coef);
    x(i+1) = pol(coef.*(2*rem(1:d,2)-1),d-1);
    if dp > 0
     coef = p{1:2:dp}; d = length(coef);
     y(i+1) = pol(coef.*(2*rem(1:d,2)-1),d-1);
    end;
  end;

  % compute numerator and denominator of transfer function
  % to be minimized to obtain stability radius
  
  % compute pseudo-inverse of polynomial matrix coef
  % first minor
  a1 = pol(0);
  for i = 1:n,
    a1 = a1 + x(i+1)^2 + y(i+1)^2;
  end;
  a1 = a1^2;
  % second minor
  sumx = pol(0); sumy = pol(0); sumxy = pol(0);
  for i = 1:n,
    sumx = sumx + x(i+1)^2;
    sumy = sumy + y(i+1)^2;
    sumxy = sumxy + x(i+1)*y(i+1);
  end;
  a2 = sumx*sumy - sumxy^2;

  % first denominator
  b1 = pol(0);
  for i = 1:n,
    b1 = b1+ (x(i+1)*x(1)+y(1)*y(i+1))^2;
  end;

  % second denominator
  b2 = pol(0);
  for i = 1:n,
    b2 = b2 + (x(i+1)*y(1)-x(1)*y(i+1))^2;
  end;
     
  d = 0;
  
  if isinf(deg(a2)), % rank = 1
   if isinf(deg(a1)), % rank = 0
     d = inf;
   else
     % rank = 1
     a = a1; b = b1;
     % compute numerator of derivative
     c = a*polyder(b) - b*polyder(a);
     if deg(c) < 1 % constant function
      d = sqrt(b{0}/a{0});
     else
      % extract positive real roots
      r = roots(c);
      r = real(r((abs(imag(r))<tol)&(real(r)>tol)));
      % evaluate transfer function at positive real roots
      % and extract minimum value
      d = inf;
      for i = 1:length(r),
	d = min(d,sqrt(polyval(b,r(i))/polyval(a,r(i))));
      end;
     end;
   end;
  else
    % rank = 2
    a = a2; b = b2;
    % compute numerator of derivative
    c = a*polyder(b) - b*polyder(a);
    if deg(c) < 1 % constant function
      d = sqrt(b{0}/a{0});
    else
      % extract positive real roots
      r = roots(c);
      r = real(r((abs(imag(r))<tol)&(real(r)>tol)));
      % evaluate transfer function at positive real roots
      % and extract minimum value
      d = inf;
      for i = 1:length(r),
	d = min(d,sqrt(polyval(b,r(i))/polyval(a,r(i))));
      end;
    end;
    % detect discontinuities
    r = roots(a);
    r = real(r((abs(imag(r))<tol)&(real(r)>tol)));
    for i = 1:length(r),
      if abs(polyval(b,r)) < tol,
	d = min(d,sqrt(polyval(b1,r(i))/polyval(a1,r(i))));
      end;
    end;
  end;
  
  % special case w = 0
  n0 = norm(pq{0});
  if n0 > 0,
    d0 = abs(p0{0})/n0;
  else
    d0 = inf;
  end;

  % special case w = inf
  dp = max(deg(p0),deg(pq));
  ninf = norm(pq{dp});
  if ninf > 0,
    dinf = abs(p0{dp})/ninf;
  else
    dinf = inf;
  end;

end;

if verbose,
 disp(['ELLIANA: stability margin to the origin = ' num2str(d0)]);
 disp(['ELLIANA: stability margin to the imaginary axis = ' num2str(d)]);
 disp(['ELLIANA: stability margin to infinity = ' num2str(dinf)]);
end;

r = min([d0 dinf d]);

