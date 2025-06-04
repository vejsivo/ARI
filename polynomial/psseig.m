function L = psseig(F, G, P, tol)
%PSSEIG  Polynomial approach to eigenstructure assignment
%         for state-space system
%
% Given a linear system
%   .
%   x = F*x + G*u
% 
% where F is an NxN constant matrix and G is an NxM constant matrix, and
% a set of polynomials P = {P{1}, P{2}, .., P{R}} where R <= M, the command
%
%   L = PSSEIG(F, G, P)
%
% returns, when possible, a constant matrix L such that the closed-loop
% matrix of the controlled system
%   .
%   x = (F-G*L)*x
%
% has invariant polynomials Q{1}, Q{2}, .., Q{N} where
%
%   Q{1} = P{1}*Q{2}, Q{2} = P{2}*Q{3}, ..,
%   Q{R} = P{R}, Q{R+1} = .. = Q{N} = 1.
%
% Such a matrix exists if and only if the fundamental degree inequality
%
%  DEG(Q{1}) + DEG(Q{2}) + .. + DEG(Q{K}) >= c1 + c2 + .. + cK
%
% holds for all K = 1, 2, .., R, where c1 >= c2 >= .. >= cR are
% the controllability indices of pair (F,G). Moreover, equality must
% hold for K = R. If input polynomials P do not satisfy these conditions,
% then the macro issues an error message.
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.

%    Author: D. Henrion, September 12, 1999.
%    Copyright 1999 by Polyx, Ltd.
%    Modified by J. Jezek, Aug 2001, arg checking

global PGLOBAL;
eval('PGLOBAL.VERBOSE;', 'painit;');

verbose = strcmp(PGLOBAL.VERBOSE, 'yes');

% ---------------------
% Parse input arguments
% ---------------------

if nargin < 3,
   error('Not enough input arguments.');
end;

if nargin < 4 | isempty(tol),
   tol = PGLOBAL.ZEROING;
else
   if ~isnumeric(tol) | length(tol)~=1 | ~isreal(tol) | ...
         tol<0 | tol>1,
      error('Invalid tolerance.');
   end;
end;

if ~isnumeric(F) | ~isnumeric(G) | ndims(F)>2 | ndims(G)>2,
   error('Invalid 1st or 2nd argument.');
end;
n = size(F, 1); m = size(G, 2);
if (size(F, 2) ~= n) | (size(G, 1) ~= n),
   error('Matrices of inconsistent dimensions.');
end;

if ~isa(P, 'cell'),
   if max(size(P)) > 1,
      error('Invalid 3rd argument.');
   else
      P = {P};
   end;
end;

P = P(:); r = length(P);

if (r > n) | (r > m),
   error('Too many invariant polynomials in 3rd argument.');
end;

% ---------------------------------
% Build monic invariant polynomials
% ---------------------------------

if verbose,
 disp('PSSEIG: Build invariant polynomials.');
end;

Q = cell(1, m);
var = ''; h = []; invar = 0; inh = 0;

if m > 0,
   for i = m:-1:1,
      if i > r,
         Q{i} = pol(1);
      else,
         eval('PP = pol(P{i});', ...
            'error(''Some entry in 3rd argument cannot be converted to polynomial.'');');
         if any(size(PP)~=[1 1]),
            error('Invalid some entry in 3rd argument; must be scalar polynomial.');
         end;
         if isinf(deg(PP)),
            error('Invalid some entry in 3rd argument; must be nonzero polynomial.');
         end;
         if isempty(var),
            var = PP.v;
         elseif ~isempty(PP.v) & ~strcmp(var,PP.v),
            if length([var,PP.v])>2,
               error('Inconsistent variables in 3rd argument.');
            end;
            PP.v = var; invar = 1;
         end;
         if isempty(h) | ~isfinite(h),
            h = PP.h;
         elseif ~isempty(PP.h) & isfinite(PP.h) & h~=PP.h,
            PP.h = h; inh = 1;
         end;
         if i == m,
            Q{i} = PP / lcoef(PP);
         else
            Q{i} = Q{i+1} * PP / lcoef(PP);
         end;
      end;
   end;
end;   

if invar,
   warning('Inconsistent variables in 3rd argument.');
elseif inh,
   warning('Inconsistent sampling periods in 3rd argument.');
end;

% ------------------------
% Computation of right MFD
% ------------------------

if verbose,
 disp('PSSEIG: Compute MFD');
end;

[BR, AR] = ss2rmf(F, G, eye(n), zeros(n,m), tol);
pprop(AR, var); pprop(BR, var);

% sorted controllability indices c(1) >= c(2) >= ..

c = -sort(-deg(AR, 'col'));

% check the fundamental degree conditions

if verbose,
 disp('PSSEIG: Check degree conditions.');
end;

if m > 0,
   for i = 1:m,
      sum1 = 0; sum2 = 0;
      for j = 1:i,
         sum1 = sum1 + deg(Q{j});
         sum2 = sum2 + c(j);
      end;
      if (sum1 < sum2) | ((i == m) & (sum1 ~= sum2)),
         error('Invariant polynomials do not satisfy the fundamental degree conditions.');
      end;
   end;
end;   

% --------------------------------------------
% Build an MxM polynomial matrix CR
% whose invariant polynomials are Q{1},..,Q{R}
% and whose column degrees are c(1),..,c(N)
% --------------------------------------------

if verbose,
 disp('PSSEIG: Build right hand-side polynomial matrix.');
end;

CR = pol(eye(m), var);

if r > 0,
   for i = 1:r,
      CR(i,i) = Q{i};
   end;
end;

c2 = deg(CR, 'col');

stop = 0;

while ~stop,

 stop = 1;  i = 1;
 while (i <= m) & stop,
  stop = c2(i) >= c(i);
  i = i + 1;
 end;
 i = i - 1;

 if ~stop,

  stop = 1; j = 1;
  while (j <= m) & stop,
   stop = c2(j) <= c(j);
   j = j + 1;
  end;
  j = j - 1;

  if stop,
   error('Invalid column degrees.');
  end;

  % switch columns

  CR(j,:) = CR(j,:) + pol([0 1], 1, var) * CR(i,:);
  alpha = CR{c2(j)}(j,j); k = c2(j) - c2(i) - 1;
  CR(:,j) = CR(:,j) - alpha * pol([zeros(1,k) 1], k, var) * CR(:,i);

  CR = pzer(lcoef(CR, 'col') \ CR, tol);
  c2 = deg(CR, 'col');

 end;

end;

% ---------------------------------------------------
% Solve the polynomial matrix equation XL*AR+YL*BR=CR
% for constant matrices XL, YL
% ---------------------------------------------------

if verbose,
 disp('PSSEIG: Solve Diophantine equation.');
end;

[XL, YL] = xaybc(AR, BR, CR, 0, tol);

% ----------------------------
% Build feedback gain matrix L
% ----------------------------

L = XL\YL;

%end .. pseig
