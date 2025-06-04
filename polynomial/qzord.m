function [AA,BB,Q,Z] = qzord(A,B,arg3,arg4);
%QZORD  Ordered qz transformation
%
% The command
%    [AA,BB,Q,Z] = QZORD(A,B[,METHOD][,TOL])
% produces the ordered qz transformation of the matrices A and B.
% AA and BB are upper triangular matrices and Q and Z nonsingular
% matrices such that Q*A*Z = AA and Q*B*Z = BB. The generalized 
% eigenvalues of the matrix pair (A,B) are the ratios of the 
% diagonal entries of AA and BB.
%
% The arguments METHOD and TOL are optional.
% - If METHOD = 'partial' then the generalized eigenvalues are ordered 
%   such that the infinite generalized eigenvalues come last. This
%   is the default method.
% - If METHOD = 'full' then the generalized eigenvalues are ordered 
%   according to increasing real parts with the infinite generalized 
%   eigenvalues last.
% The input argument TOL sets the tolerance in determining whether a 
% generalized value is infinite.

% Author: Huibert Kwakernaak, 1997
% Modified June, 2000
%          August 2001, J.Jezek,  arg chccking
%          September 2009, M. Sebek, minor fix
% Copyright 1998, 2000 by PolyX Ltd.    

% Initialization

if nargin<2,
   error('Not enough input aguments.');
end;
if ~isnumeric(A) | ndims(A)>2,
   error('Invalid 1st argument; must be a matrix.');
end;
if ~isnumeric(A) | ndims(A)>2,
   error('Invalid 2nd argument; must be a matrix.');
end;
[rA,cA] = size(A); [rB,cB] = size(B);
if rA ~= cA | rB ~= cB | rA ~= rB
   error('Matrices of inconsistent dimensions.')
end

method = 'partial'; epp = eps;
if nargin >= 3
   if isstr(arg3)
      method = arg3;
   elseif isnumeric(arg3)
      epp = arg3;
   else
      error('Invalid 3rd argument.');
   end
end
if nargin == 4
   if isstr(arg4)   
      method = arg4;
   elseif isnumeric(arg4)  %corrected by MS, it was 'numeric'
      epp = arg4;
   else
      error('Invalid 4th argument.');
   end;
end
if ~(strcmp(method,'partial') | strcmp(method,'full'))
   error('Invalid command option.')
end
if length(epp)~=1 | ~isreal(epp) | epp<0 | epp>1,
   error('Invalid tolerance.');
end;
n = rA;

[AA,BB,Q,Z] = qz(A,B);
if length(A) == 0, return, end

% Order the diagonal entries of AA/BB 

a = diag(AA); b = diag(BB);
for i = 1:n
   if abs(b(i)) < epp*norm(BB)*max(size(BB))^2*abs(a(i))
      d(i) = Inf;
   else
      if strcmp(method,'partial')
         d(i) = 0;
      else
         d(i) = real(a(i)/b(i));
      end
   end
end
d = d';  

for i = 1 : n*n/4
% [a b d],keyboard
   if sort(d) == d
      break               % done ordering
   end
   for k = 1:n-1
      a11 = AA(k,k); a12 = AA(k,k+1); a22 = AA(k+1,k+1);
      b11 = BB(k,k); b12 = BB(k,k+1); b22 = BB(k+1,k+1);
      if d(k) > d(k+1)
         [c1,s1] = cgivens1(a11*b12-a12*b11,a11*b22-a22*b11,epp);
         [c2,s2] = cgivens1(a12*b22-a22*b12,a11*b22-a22*b11,epp);
         AA(:,k:k+1) = AA(:,k:k+1)*[c2 s2;-s2' c2];
         AA(k:k+1,:) = [c1 s1;-s1' c1]*AA(k:k+1,:);
         BB(:,k:k+1) = BB(:,k:k+1)*[c2 s2;-s2' c2];
         BB(k:k+1,:) = [c1 s1;-s1' c1]*BB(k:k+1,:);
         Z(:,k:k+1) = Z(:,k:k+1)*[c2 s2;-s2' c2];
         Q(k:k+1,:) = [c1 s1;-s1' c1]*Q(k:k+1,:);
         tt = a(k); a(k) = a(k+1); a(k+1) = tt;
         tt = b(k); b(k) = b(k+1); b(k+1) = tt;
         tt = d(k); d(k) = d(k+1); d(k+1) = tt;
      end
   end
end

%end .. qzord
