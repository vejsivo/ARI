function in = inertia(A,tol)
%INERTIA  Inertia of a polynomial matrix
%
% The command
%    IN = INERTIA(A,TOL) 
% returns the inertia of the polynomial matrix A. The inertia is 
% the triple [ns n0 nu], with ns the number of stable roots of A, 
% n0 the number of marginally stable roots, and nu the number of 
% unstable roots.
%
% The definition of stability depends on the variable of the 
% polynomial matrix. The root r of a polynomial matrix is 
% considered as stable
%  . with variable 's' or 'p':    if real(r) < 0
%  . with variable 'z^-1' or 'd': if abs(r) > 1
%  . with variable 'z' or 'q':    if abs(r) < 1
%
% The optional input parameter TOL is a tolerance which is used to 
% determine if certain coefficients that arise in the construction 
% of the Routh array are zero. Its default value is the global 
% zeroing tolerance.
%
% The algorithm and its code are from
% G. Meinsma, "Elementary proof of the Routh-Hurwitz test."
% Systems & Control Letters 25 (1995) 237-242.

%    Author: H. Kwakenaak, August, 1998. Modified January, 1999.
%    Copyright 1998 by Polyx, Ltd.


% Intialization and checks

global PGLOBAL;
eval('PGLOBAL.ZEROING;', 'painit;');

if nargin == 0,
   error('Not enough input arguments.');
end;
if nargin == 1 | isempty(tol),
   tol = PGLOBAL.ZEROING;
else
   if ~isa(tol,'double') | length(tol)~=1 | ...
         ~isreal(tol) | tol<0 | tol>1,
      error('Invalid tolerance.');
   end;
end;

eval('A = pol(A);', 'error(peel(lasterr));');
[ra, ca] = size(A);
type = A.var;
if ra ~= ca,
 error('Matrix is not square.'); 
end;
if isempty(A)
   in = []; return
end


% Compute the determinant of A if necessary

if (ra > 1) | (ca > 1),
   a = det(A);
else
   a = A;
end;


% Convert to the continuous-time case

if strcmp(type,'z^-1') | strcmp(type,'d')   % conjugate
 a = a'; pprop(a,'z'); type = 'z';
end;
if strcmp(type,'z')  % bilinear transformation
 aa = 0; 
 for i = 0:deg(a)
  aa = aa+a{i}*(1-z)^(deg(a)-i)*(1+z)^i;
 end
 pprop(aa,'s'); a = aa;
end


% Zero any small coefficients

a = pzer(a,tol*norm(a));


% Constant matrix case

dega = deg(a);
if dega <= 0
   in = [0 0 0]; return
end


% Compute the inertia

a = a{dega:-1:0};
in = [0 0 0];
wehavehadcase2 = 0;
for n = dega:-1:1   % Reduce the degree to 1
    k = find(abs(a(2:2:n+1)) > tol);
    if isempty(k)   % Case 2: Differentiate
       a(2:2:n+1) = a(1:2:n).*(n:-2:1);
       wehavehadcase2 = 1;
    elseif k(1) > 1 % Case 1: Add polynomial
       ind = 0:2:(n+1-2*k(1));
       f = (-1)^k(1);
       a(ind+2) = a(ind+2)-f*a(ind+2*k(1));
    end
    eta = a(1)/a(2);
    if wehavehadcase2
         in = in+[(eta<0) 0 (eta<0)];
    else in = in+[(eta>0) 0 (eta<0)];
    end
    a(1:2:n) = a(1:2:n)-eta*a(2:2:n+1);
    a(1) = [];       % Reduce degree to n-1
end
in = in+[0 dega-sum(in) 0];

%end .. inertia

