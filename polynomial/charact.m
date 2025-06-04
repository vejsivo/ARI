function [V,N,G,M] = charact(P, Z, tol)
%CHARACT  Characteristic vectors of a polynomial matrix
%
% Given a nonsingular polynomial matrix P and a complex number Z,
% the command
%    [V,N,G,M] = CHARACT(P, Z) 
% returns a cell array V of generalized right characteristic vectors of P
% associated with the zero Z. The zero Z has algebraic multiplicity N and
% geometric multiplicity G.
%    
% . If Z has algebraic multiplicity N = M(1) + .. + M(G) then there are N
%   nonzero generalized characteristic vectors V{1,1:M(1)},..,V{G,1:M(G)}.
%   The vectors are ordered in G chains V{I,:}, each of length M(I) and
%   such that M(1) >= .. >= M(G). The first vectors in each chain, that is,
%   the vectors V{1:G,1}, are linearly independent.
%
%   Let P[K] denote the Kth derivate of P evaluated at Z.
%   Then the characteristic vectors in chain number I satisfy
%
%      0 = P[0]*V{I,1}
%      0 = P[0]*V{I,2}  + P[1]*V{I,1}
%      ...
%      0 = P[0]*V{I,M(I)} + P[1]*V{I,M(I)-1} +..+ P[M(I)-1]*V{I,1}/(M(I)-1)!
%
% . If Z is not a zero of P then V is empty and N = G = 0.
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.

%   Author: D. Henrion, November 9, 1998
%   Last modified by D. Henrion, December 16, 1999.
%   Updated to 3.0 by D. Henrion, August 30, 2000.
%   Copyright 1998-2000 by Polyx, Ltd.
%

%   The macro is based on succesive extractions by CEF and NULLREF
%   of right null spaces of constant Toeplitz matrices
%
%   T = 
%
%     [ P[0]  P[1]  P[2]/2!
%       0     P[0]  P[1]
%       0     0     P[0]
%                              ..]

global PGLOBAL;
eval('PGLOBAL.VERBOSE;', 'painit;');
verbose = strcmp(PGLOBAL.VERBOSE, 'yes');

if nargin<2,
   error('Not enogh input arguments.');
end;
eval('P = pol(P);', 'error(peel(lasterr));');
[n, m] = size(P);
degP = max(deg(P),0);

if n ~= m,
 error('Matrix is not square.');
end;

if ~isa(Z, 'double') | (length(Z) > 1),
   error('Invalid 2nd argument.');
end;
if nargin < 3 | isempty(tol),
   tol = PGLOBAL.ZEROING;
else
   if ~isa(tol,'double') | length(tol)>1 | ...
         ~isreal(tol) | tol<0 | tol>1,
      error('Invalid tolerance.');
   end;
end;

if rank(P, tol) < n,
 error('Matrix is singular.');
end;

% Successive extractions of null-spaces of Toeplitz matrix T

if verbose,
  disp('CHARACT: Toeplitz matrix null-space extraction.');
end;

T = polyval(P, Z);
T(abs(T) < tol) = 0;

stop = 0;
order = 0; chains = [];

if verbose,
 fprintf('CHARACT: Toeplitz order = ');
end;  

while ~stop,
  
  order = order + 1;

  if verbose,
   fprintf('%d ', order);
  end;

  % Detect linearly dependent columns in Toeplitz matrix T
  % by reduction into row echelon form

  [K,icol] = cef(T',tol);

  rankK = length(icol);
  depcol = ones(1, order*n); depcol(icol) = 0;
  newdepcol = depcol(1+(order-1)*n:order*n);
  stop = sum(newdepcol) == 0; % no new dependent columns ?

  % Detect new independent columns:
  % They correspond to new chains of characteristic vectors

  if order > 1,

   olddepcol = depcol(1+(order-2)*n:(order-1)*n);
   newdepcol = xor(newdepcol, olddepcol);

   for i = 1:n,
    if newdepcol(i),
     % store column index in null-space
     chains = [chains [order-1; sum(depcol(1:i+(order-2)*n))]];
    end;
   end;

  end;

  if ~stop,

   if order <= degP*n,

    % Append to Toeplitz matrix T next derivative of P evaluated at Z

    P = polyder(P) / order;
    newT = polyval(P, Z);
    newT(abs(newT) < tol) = 0;
    T = [[T(1:n,:) newT]; [zeros(order*n,n) T]];

   else
   
    % Stop process due to numerical problems
    % Return empty output argument
    stop = 1;

   end;

  end;

end;

if verbose,
 fprintf('\n');
end;

% Geometric multiplicity: number of chains of char. vectors

G = size(chains, 2);
M = zeros(G, 1);
N = 0;

V = {};

if G > 0,

  if verbose,
   disp('CHARACT: Compute characteristic vectors.');
  end;

  % Compute char. vectors by null-space extraction

  K = K(1:end-n, 1:end-n);
  icol = icol(1:end-n);

  K = nullref(K, icol)';

  for i = 1:G,

   k = chains(1, i); % length of chain
   l = chains(2, i); % column index in null-space

   % Store char. vectors

   % Make sure there is no zero vector in the chain

   j = 0;
   while (j < k),

    j = j + 1;

    vector = K(1+(k-j)*n:(k-j+1)*n, l);

    if norm(vector) < tol, % zero vector ?

     if verbose,
      disp('CHARACT: Found a zero vector in the chain.');
     end;

     % Use next non-zero column in null-space
     l2 = l-1; vector = K(1+(k-j)*n:(k-j+1)*n, l2);
     while (norm(vector) < tol) & (l2 > 0),
      vector = K(1+(k-j)*n:(k-j+1)*n, l2);
      l2 = l2 - 1;
     end;

     if l2 == 0,

      error('CHARACT: Invalid zero characteristic vector.');

     else

      % linear combination of the columns to make vector non-zero

      K(:, l) = K(:, l) + K(:, l2);

     end; 

    end;

   end;

   % Extract characteristic vectors in the chain

   for j = 1:k,

    N = N + 1;
    V{G-i+1, j} = K(1+(k-j)*n:(k-j+1)*n, l);

   end;

   M(G-i+1) = k;

  end;

end;

%end .. charact
