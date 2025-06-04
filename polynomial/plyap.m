function [X,Y] = plyap(A,B,C,tol)
%PLYAP  Solution of a two-sided equation with matrix pencil coefficients
%
% The command
%    [X,Y] = plyap(A,B,C)
%    [X,Y] = plyap(A,B,C,tol)
% solves the "pencil Lyapunov equation"
%    A*X + Y*B = C
% with A,B and C polynomial matrices of degree 1, A,B being square,
% for the constant matrices X and Y. 
%
% If A and B have no common zeros then a unique solution exists. If 
% this condition is not satisfied then the solution may be inaccurate.
%
% The optional parameter tol is a tolerance that is used in testing 
% whether any of the pencils A or B is upper or lower triangular.
%
% If the verbose property is set to 'yes' then the macro reports the
% relative error.

% Author: Huibert Kwakernaak, 1996-1998
% Copyright 1998 by PolyX, Ltd
% Modified by Jan Jezek, Aug 2001, Jan 2002, arg checking


% Initialize and convert to the equation
%  (s*E1-A1)*X + Y*(s*E2-A2) = s*E3-A3

global PGLOBAL;
eval('PGLOBAL.ZEROING;', 'painit;');

if nargin<3,
   error('Not enough input arguments.');
end;
eval('A = pol(A); B = pol(B); C = pol(C);', ...
   'error(peel(lasterr));');

if size(A,1)~=size(A,2) | A.deg>1
   error('1st argument is not a square pencil.'); 
else
   n1 = size(A,1); 
   A1 = -A{0}; 
   if A.deg == 1
      E1 = A{1};
   else
      E1 = zeros(size(A));
   end
end

if size(B,1)~=size(B,2) | B.deg>1
   error('2nd argument is not a square pencil.'); 
else
   n2 = size(B,2); 
   A2 = -B{0}; 
   if B.deg == 1
      E2 = B{1};
   else
      E2 = zeros(size(B));
   end
end

if  size(C,1)~=size(A,1) | size(C,2)~=size(B,2) | C.deg>1
   error('3rd argument is not a pencil of correct dimensions.')
else
   A3 = -C{0}; 
   if C.deg == 1
      E3 = C{1};
   else
      E3 = zeros(size(C));
   end
end

[tv,v,A,B,C] = testvp3(A,B,C);
if tv==2,
   error('Inconsistent variables.');
elseif tv==0,
   warning('Inconsistent variables.');
end;

if nargin == 3 | isempty(tol),
   tol = sqrt(eps);
else
   if ~isa(tol,'double') | length(tol)~=1 | ...
         ~isreal(tol) | tol<0 | tol>1,
      error('Invalid tolerance.');
   end;
end

EE1 = E1; AA1 = A1; EE2 = E2; AA2 = A2; EE3 = E3; AA3 = A3;


% Check if the problem has any complex inputs

if any(any(imag(E1))) | any(any(imag(A1))) ...
   | any(any(imag(E2))) | any(any(imag(A2))) ...
   | any(any(imag(E3))) | any(any(imag(A3)))
   real_flg = 0;
else
   real_flg = 1;
end


% Test if s*E1-A1 or s*E2-A2 is upper or lower triangular

case1 = 0;  % Neither pencil is upper or lower triangular

e = zeros(n1); a = e;
for i = 1:n1, for j = 1:i-1
   e(i,j) = E1(i,j); a(i,j) = A1(i,j);
end, end
if norm(e,1) < tol*norm(E1) & norm(a,1) < tol*norm(A1)
   case1 = 1;  % s*E1-A1 is upper triangular
end

e = zeros(n1); a = e;
for i = 1:n1, for j = i+1:n1
   e(i,j) = E1(i,j); a(i,j) = A1(i,j);
end, end   
if norm(e,1) < tol*norm(E1) & norm(a,1) < tol*norm(A1)
   case1 = 2;  % s*E1-A1 is lower triangular
end 

e = zeros(n1); a = e;
for i = 1:n2, for j = 1:i-1
   e(i,j) = E2(i,j); a(i,j) = A2(i,j);
end, end
if norm(e,1) < tol*norm(E2) & norm(a,1) < tol*norm(A2)
   case1 = 3;  % s*E2-A2 is upper triangular
end

e = zeros(n2); a = e;
for i = 1:n2, for j = i+1:n2
   e(i,j) = E2(i,j); a(i,j) = A2(i,j);
end, end   
if norm(e,1) < tol*norm(E1) & norm(a,1) < tol*norm(A1)
   case1 = 4;  % s*E2-A2 is lower triangular
end 


% Arrange s*E2-A2 to be upper triangular

if case1 == 0
   [E2,A2,Q,Z] = qz(E2,A2); E3 = E3*Z; A3 = A3*Z;
elseif case1 == 1          % s*E-A1 is UT
   e = E1; a = A1;
   E1 = E2'; A1 = A2'; E3 = fliplr(E3'); A3 = fliplr(A3');
   E2 = flipud(fliplr(e')); A2 = flipud(fliplr(a')); 
elseif case1 == 2          % s*E1-A1 is LT
   e = E1; a = A1; 
   E1 = E2'; A1 = A2'; E3 = E3'; A3 = A3';
   E2 = e'; A2 = a';
elseif case1 == 4          % s*E2-A2 is LT
   E2 = flipud(fliplr(E2)); A2 = flipud(fliplr(A2)); 
   E3 = fliplr(E3); A3 = fliplr(A3);
end
n1 = length(E1); n2 = length(E2);


% Solve for X and Y column by column from left to right

X = zeros(n1,n2); Y = X;

for j = 1:n2
   e = E2(j,j); a = A2(j,j);
   X(:,j) = (a*E1-e*A1)\(a*E3(:,j)-e*A3(:,j)-Y*(a*E2(:,j)-e*A2(:,j)));
   Y(:,j) = e'*E3(:,j)+a'*A3(:,j)-(e'*E1+a'*A1)*X(:,j) ...
            -Y*(e'*E2(:,j)+a'*A2(:,j));
   Y(:,j) = Y(:,j)/(abs(e)^2+abs(a)^2);
end


% Reconstruct the solution

if case1 == 0
   X = X/Z; Y = Y*Q;
elseif case1 == 1
   X = fliplr(X); Y = fliplr(Y);
   x = X; X = Y'; Y = x';
elseif case1 == 2
   x = X; X = Y'; Y = x';
elseif case1 == 4
   X = fliplr(X); Y = fliplr(Y);
end


% Make the result real

if real_flg == 1
   X = real(X); Y = real(Y);
end


% Verbose level

verbose = strcmp(PGLOBAL.VERBOSE, 'yes');


% Checks

Eps0 = norm(A*X+Y*B-C,'blk',1);
Nrm = max([norm(A,'blk',1),norm(B,'blk',1),norm(C,'blk',1)]);
Eps = Eps0/Nrm;
if verbose
   disp(sprintf('plyap: Relative residue %g',Eps));
elseif Eps > 1e-6
   disp(sprintf('plyap warning: Relative residue %g',Eps));
end
if any(any(isnan(X))) | any(any(isinf(X))) | ...
      any(any(isnan(Y))) | any(any(isinf(Y)))
   disp(sprintf...
      ('plyap warning: Solution has Not-a-Number or Infinite entries.'));  
end

%end .. plyap
