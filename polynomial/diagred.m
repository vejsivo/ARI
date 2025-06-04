function [Zred,U,UI] = diagred(Z,tol)
%DIAGRED  Diagonal reduction
%
% If Z is a continuous-time para-Hermitian polynomial matrix then
%    [ZR,U,UI] = DIAGRED(Z)
% produces a diagonally reduced para-Hermitian matrix ZR and a 
% unimodular matrix U such that
%    ZR = U'ZU
% UI is the inverse of U.
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.

%    Author: D. Henrion, January 27, 1999.
%    Modified by D. Henrion, February 16, 2000.
%    Updated to 3.0 by D. Henrion, August 30, 2000.
%                   by J. Jezek, 10-Feb-2003  arg check
%    Copyright 1999-2000 by Polyx, Ltd.

global PGLOBAL;

if nargin < 1,
   error('Not enough input arguments.');
end;
if isa(Z,'tsp'),
   error('Only continuous-time.');
end;clear fun
eval('Z = pol(Z);', 'error(peel(lasterr));');

[n, cZ] = size(Z);
typeZ = Z.var;

if ~strcmp(typeZ, 's') & ~strcmp(typeZ, 'p') & ~isempty(typeZ),
 error('Only continuous-time.');
end;

if nargin == 1 | isempty(tol),
   tol = PGLOBAL.ZEROING;
else   
   if ~isa(tol,'double') | length(tol)~=1 | ~isreal(tol) | ...
         tol<0 | tol>1,
      error('Invalid tolerance.');
   end;
end;

me = exp(log10(1+norm(Z))); % relative tolerance
tolzero = tol * me;

if n ~= cZ,
 error('Matrix is not square.');
end;

if issingular(Z, tol),
 error('Matrix is singular.');
end;

U = pol(eye(n), 0); UI = U;
Zred = Z; reduced = 0;

% half diagonal degrees
D = 0;
eval('D = deg(Zred,''dia'');', 'error(peel(lasterr));');

while ~reduced,

 % diagonal leading coefficient corresponding
 % to half diagonal degrees D(1),..,D(n)
 degZ = max(0, deg(Zred));
 ZL = zeros(n);
 for i = 1:n, for j = 1:n,
  if D(i)+D(j) <= degZ,
   ZR = Zred{D(i)+D(j)};
   ZL(i, j) = (-1)^D(i) * ZR(i,j);
  end;
 end; end;

 reduced = (min(svd(ZL)) > tolzero) | ~any(D);

 if ~reduced, % if matrix is not diagonally reduced

  % extract right null-space of diagonal leading coefficient
  [void,S,V] = svd(ZL,0);
  if n > 1, s = diag(S);
  elseif n == 1, s = S(1);
  else s = 0;
  end;
  r = sum(s > tolzero);  
  N = V(:,r+1:n);

  deflation = 0;

  reduced = 1;
  i = 1; 
  while (i <= size(N,2)),
   e = N(:, i);
   hideg = []; maxi = -1;
   for i = 1:n,
    if abs(e(i)) > tolzero, % active index set
     if D(i) >= maxi,
      if D(i) > maxi,
       maxi = D(i); hideg = i;
      else
       hideg = [hideg i]; % highest degree active index set
      end;
     end;
    end;
   end;

   if isempty(hideg),
     error('The diagonal reduction failed. Try to modify the tolerance.');
   end;

   k = 0; maxi = -1;
   for i = hideg,
    if abs(e(i)) > maxi,
     maxi = abs(e(i)); k = i; % index corresponding to maximal element
    end;
   end;

   % deflation
   if D(k) > 0,
    reduced = 0;
    a = e/e(k);

    % build matrix U(s) and its inverse UI(s)
    V2 = pol(eye(n), 0); V2I = V2;
    for i = 1:n,
     if i ~= k,
      d = D(k) - D(i);
      if d >= 0,
       V2(i,k) = pol([zeros(1,d) a(i)], d, typeZ);
       V2I(i,k) = pol([zeros(1,d) -a(i)], d, typeZ);
      end;
     end;
    end;
    Zred = pzer(V2'*Zred*V2, tolzero); 
    U = pzer(U*V2, tolzero);
    UI = pzer(V2I*UI, tolzero);

    % new half diagonal degree
    degD = max(0, deg(Zred(k,k)));
    if rem(degD, 2),
     error('The diagonal reduction failed. Try to modify the tolerance.');
    else
     D(k) = degD / 2;
    end;

   end; % if D(k)

   i = i + 1;
  end;
  
 end; % if reduced

end; % while

%end .. diagred



