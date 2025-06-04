function s = isstable(A,tol)
%ISSTABLE  Test if polynomial matrix is stable
%
% ISSTABLE(A) returns 1 if the square polynomial matrix A is
% stable and 0 otherwise.
%
% The definition of stability depends on the variable of the
% polynomial matrix. A polynomial matrix with roots Zi is 
% considered to be stable
%    . with variable 's' or 'p':    if REAL(Zi) < 0 for all i
%    . with variable'z^-1' or 'd':  if ABS(Zi) > 1 for all i
%    . with variable 'z' or 'q':    if ABS(Zi) < 1 for all i
% In all other cases the matrix is considered to be unstable.
%
% If the matrix is empty or constant or unimodular then it is
% considered to be stable. If the matrix is singular then it is
% considered to be unstable.
%
% An absolute tolerance TOL may be specified as an additional input
% argument. It is used in testing whether an element of the Routh array
% is nonpositive. Its default value is the global zeroing tolerance.

%    Authors: D. Henrion, J. Jezek, January 11, 1999.
%    Last modified by D. Henrion, January 8, 2003.
%    Copyright 1999-2003 by Polyx, Ltd.

global PGLOBAL;

ni = nargin;
if ni<1,
   error('Not enough input arguments.');
end;
if ni == 1 | isempty(tol),
   tol = PGLOBAL.ZEROING;
else
   if ~isa(tol,'double') | length(tol)~=1 | ...
         ~isreal(tol) | tol<0 | tol>1,
      error('Invalid tolerance.');
   end;
end;

eval('A = pol(A);', 'error(peel(lasterr));');

[ra, ca] = size(A);
if ra ~= ca,
 error('Matrix is not square.');
end;
if ra == 0,
   % empty : matrix is considered as stable
   s = logical(1); return;
end;

% compute the determinant of A if necessary
if (ra > 1) | (ca > 1),
   a = det(A, tol);
else
   a = A;
end;
dega = deg(a);

if dega == 0,
   % constant non-zero determinant: matrix is considered as stable
   s = logical(1); return;
elseif dega < 0,
   % zero determinant: matrix is considered as unstable 
   s = logical(0); return;
end;

s = logical(0);
type = a.v;

if strcmp(type, 's') | strcmp(type, 'p'),

  % ***************
  % CONTINUOUS-TIME
  % ***************

  % Scale polynomial (does not affect stability)
  a = scale(a);
  
  % Complex Routh array on det(A)
  % unit constant coefficient
  if abs(a{0}) > tol
   a{:} = a{:} / a{0};
  else
   return;
  end;

  % symmetric and skew-symmetric parts
  A = zeros(2,dega+1);
  for i = 0:dega,
   c = a{i};
   if rem(i, 2),
    A(1, i+1) = imag(c); A(2, i+1) = real(c);
   else
    A(1, i+1) = real(c); A(2, i+1) = imag(c);
   end;
  end;

  l = dega+1;
  while l > 1,

   B = zeros(1,l);
   newA = zeros(2,l);

   % compute new symmetric part
   for i = 1:2:l-1,
    newA(1,i) = A(1,1)*A(2,i+1)+A(2,1)*A(1,i+1);
   end;
   for i = 2:2:l-1,
    newA(1,i) = A(1,1)*A(2,i+1)-A(2,1)*A(1,i+1);
   end;

   % compute intermediate row
   for i = 1:2:l,
    B(i) = A(1,i)*A(2,2)+A(2,i)*A(1,2);
   end;
   B(2) = 0;
   for i = 4:2:l,
    B(i) = A(1,i)*A(2,2)-A(2,i)*A(1,2);
   end;

   % compute new skew-symmetric part

   newA(2,1:(l-1)) = B(2:l) - newA(1,2:l);

   % next array
   A = newA;
   l = l - 1;

   % stopping criterion: first element in symmetric part is <= 0
   if A(1, 1) <= tol
    return;
   end;

   A = A / A(1,1);

  end; % while l > 1

  s = logical(1);

else

  % *************
  % DISCRETE-TIME
  % *************

  dega = deg(a); 
  acoef = a{:};
  if strcmp(type, 'z^-1') | strcmp(type, 'd'),
   acoef = fliplr((acoef').'); % conjugate
  end;

  % proceed by polynomial reductions:
  % discrete-time version of complex Routh array on the determinant

  while dega > 0,

   if abs(acoef(1)) - abs(acoef(dega+1)) > -tol, % check nec. cond.
    return; % polynomial is instable
   end;

   if dega >= 1,
    temp = acoef(1:dega+1);
    temp = acoef(dega+1)'*temp - acoef(1)*fliplr((temp').');
    acoef(1:dega) = temp(2:dega+1); acoef(dega+1) = 0; % drop degree
    while (dega > 0) & abs(acoef(dega+1)) < tol,
     dega = dega - 1;
    end;
    if abs(acoef(dega+1)) >= tol,
     acoef = acoef / acoef(dega+1);
    end;
   end;

  end;  

  s = logical(1); % polynomial is stable

end;

%end .. isstable
