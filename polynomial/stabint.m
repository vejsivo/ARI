function [rmin, rmax] = stabint(varargin)
%STABINT  Stability interval of a polynomial matrix
%
% Given N+1 polynomials P0, P1, .., PN with real coefficients
% such that P0 is stable and such that DEG(P0) >= MAX(DEG(P1),..,DEG(PN))
% the function
%    [Rmin, Rmax] = STABINT(P0, P1, .., PN)
% returns the largest stability interval of the uncertain polynomial
%      P(R) = P0 + P1*R + ... + PN*R^N
% so that P(R) remains stable for any real R in the interval [Rmin, Rmax].
%
% If the input arguments P0, P1, .., PN are polynomial matrices then
% the macro first transforms the matrix problem to a scalar problem
% by computing the scalar polynomial two-dimensional determinant
%     PP(R) = DET(P(R)) = PP0 + PP1*R + .. + PPN*R^N.
% Note that the macro will work properly only if the relation
%   DEG(PP0) >= MAX(DEG(PP1), .., DEG(PPN)) holds, which is not
% necessarily the case even if DEG(P0) >= MAX(DEG(P1), .., DEG(PN)).
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.

%   D. Henrion, October 30, 1998.
%   Last modified by D. Henrion, September 30, 2002.
%   Copyright 1998-2002 by Polyx, Ltd.

% The macro is based on continuous-time and discrete-time guardian maps
% for uncertain polynomials. Stability bounds are computed as polynomial
% matrix real zeros nearest to the origin. When input arguments are
% polynomial matrices, then macro DET2D is called first to transform
% the problem into a scalar polynomial problem.

global PGLOBAL;

nargin = length(varargin);

% Parse input arguments.
% Arrange input matrices into a cell array P = {P1, P2, ..}.

p = {}; n = 0; tol = [];

for i = 1:nargin,
 arg = varargin{i}; 
 ismat = 0;
 if isa(arg, 'double'),
  if (i <= 2),
   ismat = 1;
  else
   if any([1 1]-size(arg)),
    ismat = 1;
   else
    tol = arg;
   end;
  end;
 elseif isa(arg, 'pol'),
  ismat = 1;
 else
  error('Invalid input argument.');
 end;
 if ismat,
  arg = pol(arg);
  if any(any(isnan(arg))),
   error('Invalid Not-a-Number entry in the input matrices.');
  elseif any(any(isinf(arg))),
   error('Invalid infinite entry in the input matrices.');
  elseif ~isempty(arg),
   n = n + 1;
   if n == 1,
    [m,r] = size(arg);
    var = arg.var;
    d = deg(arg);
   else
    [q,r] = size(arg);
    if q ~= m,
     error('The input matrices have inconsistent dimensions.');
    end;
    if isempty(var),
     var = arg.var;
    elseif ~isempty(arg.var) & ~strcmp(var, arg.var),
     warning('Inconsistent variables.');
     var = PGLOBAL.VARIABLE;
    end;
   end;
   if m ~= r,
    error('The input matrices must be square.');
   end;
   p{n} = arg;
  end;
 end;
end;

if isempty(tol),
 tol = PGLOBAL.ZEROING;
end;

if n < 2,
 error('Two polynomial matrices or more must be specified.');
end;

if m > 1, % polynomial matrices

 % compute 2D determinant
 dc = det2d(p{:}, tol);

 % retrieve scalar polynomials 
 n = length(dc);
 p = cell(n, 1);
 for i = 1:n,
  p{i} = dc(i);
 end;

end;

% nominal polynomial

p0 = pol(p{1});
d = deg(p0);
sp = symbol(p0);

% check degrees of remaining polynomials

for i = 2:n,
 if deg(p{i}) > d,
   error('Relation on polynomial degrees does not hold. See the documentation.');
  end;
end;

switch sp,
case {'s', 'p'}, % continuous-time polynomial

 % check whether p0 is stable

 if ~isstable(p0, tol),
  error('The first polynomial must be stable.');
 end;

 % build polynomial matrix H(R) = H0 + H1 * R + .. + Hn * R^n
 % where Hi is the Hurwitz matrix of polynomial pi

 H = zeros(d, d*n);
 for i = 1:n,
  H(:, 1+(i-1)*d:i*d) = hurwitz(p{i}, d);
 end;
 H = pol(H, n-1);
 
 % roots of H(R)
 rootH = roots(H, tol);

 % roots due to possible degree drop
 drop = zeros(1,n);
 for i = 1:n,
  drop(i) = p{i}{d};
 end;
 drop = pol(drop,n-1);
 
 rootH = [rootH; roots(drop)];
 
otherwise, % discrete-time polynomial
 
 % forward of backward operator
 flag = (strcmp(sp, 'z^-1') | strcmp(sp, 'd'));

 % check whether p{0} is stable
 if ~isstable(p0),
  error('The first polynomial must be stable.');
 end;

 % detect instability at 1: build polynomial at 1
 H = zeros(1, n);
 for i = 1:n,
  if deg(p{i}) > d,
   error(['Incorrect degree of polynomial #' int2str(i) '.']);
  end;
  H(i) = polyval(p{i}, 1);
 end;
 H = pol(H, n-1);
 rootH1 = roots(H, tol);

 % detect instability at -1: build polynomial at -1
 H = zeros(1, n);
 for i = 1:n,
  H(i) = polyval(p{i}, -1);
 end;
 H = pol(H, n-1);
 rootHm1 = roots(H, tol);

 % detect instability elsewhere on the unit circle
 % with the help of the Jury matrix
 if d > 1,

  % build polynomial matrix H(R) = H0 + H1 * R + .. + Hn * R^n
  % where Hi is the Jury matrix of polynomial pi
  H = zeros(d-1, (d-1)*n);
  for i = 1:n,
   H(:, 1+(i-1)*(d-1):i*(d-1)) = jury(p{i}, d, flag);
  end;
  H = pol(H, n-1);

  % roots of H(R)
  rootHR = roots(H, tol);

 else

  rootHR = [];

 end;

 rootH = [rootH1; rootHm1; rootHR];

end;

% compute real eigenvalues nearest the origin
 
rootH = real(rootH(abs(imag(rootH)) < tol));

rmin = max(rootH(rootH < 0));
if isempty(rmin),
  rmin = -Inf;
end;

rmax = min(rootH(rootH > 0));
if isempty(rmax),
  rmax = +Inf;
end;

function S = jury(p, d, flag)
% Given a polynomial p(z) = p(0) + p(1)*z + .. + p(d)*z^d,
% compute Jury's (d-1)x(d-1) matrix 
%
% S(p) = [ p(d)  p(d-1) p(d-2) .. p(3)         p(2)-p(0)
%          0     p(d)   p(d-1) .. p(4)-p(0)    p(3)-p(1)
%                              ..
%          0     -p(0)  -p(1)  .. p(d)-p(d-4)  p(d-1)-p(d-3)
%          -p(0) -p(1)  -p(2)  ..     -p(d-3)  p(d)-p(n-2)   ]
%
% if flag == 1, coefficient p(0), p(1), .., p(d) are reversed

q = p{:};
if length(q) < d+1,
 q = [q zeros(1, d+1-length(q))];
end;

if flag,
 S = fliplr(hankel(q(d+1:-1:3))) - rot90(hankel(q(1:d-1)), 2);
else
 S = fliplr(hankel(q(3:d+1))) - rot90(hankel(q(d-1:-1:1)), 2);
end;
