function C = hermcoef(S,n)
% Given a 2x2 Hermitian matrix S and an integer N, the instruction
%
%   C = HERMCOEF(S,N)
%
% builds coefficients of the Hermite matrix of a polynomial of degree N
%
% Coefficients are stored in a matrix C such that entries of the NxN
% Hermite matrix H = [Hij] of a polynomial p(s) = p0+p1*s+p2*s^2+..
% satisfy the relationship:
%
%    [H11;H12;H22;H13;H23;H33..] = C*[p0*p0;p1*p0;p1*p1;p2*p0;p2*p1;p2*p2..]
%
% Notice that the Hermite matrix is bilinear in polynomial coefficients
%  
% Hermite matrix H is positive definite if and only if polynomial p(s)
% has its roots in the stability region described by the quadratic
% scalar inequality S(1,1)+S(1,2)*s+S(2,1)*s'+S(2,2)*s*s' < 0

% The relationship between entries Hij and coefficients pi*pj is given in
% [D. Henrion, D. Peaucelle, D. Arzelier, M. Sebek. Ellipsoidal
% Approximation of the Stability Domain of a Polynomial. Proceedings of
% the European Control Conference, pp. 384-389, Porto, Portugal, 2001]

% Author: Didier Henrion, September 17, 2001
% Last modified by Didier Henrion, November 18, 2002
% Copyright 2002 by PolyX, Ltd
  
global PGLOBAL;

eval('PGLOBAL.FORMAT;',...
     'error(''Use PINIT to initialize the Polynomial Toolbox.'');');

if nargin < 2,
 error('Invalid number of input arguments');
end;

if n < 1,
 error('Second input argument must be a positive integer');
end;

R = [eye(n) zeros(n,1); zeros(n,1) eye(n)]; % projection matrix

% Build LS A*h = B*pp with h = [H11;H21;H22;H31;H32;H33..]
% and pp = [p0*p0;p1*p0;p1*p1;p2*p0;p2*p1;p2*p2..]
% Then solve the LS with the SVD

% Build left hand-side matrix A column-wise
A = zeros((n+1)*(n+2)/2,n*(n+1)/2);
colA = 1;
for row = 1:n,
  for col = 1:row,
    H = zeros(n); H(row,col) = 1; H(col,row) = 1;
    mat = R'*kron(S,H)*R;
    A(:, colA) = trivec(mat);
    colA = colA+1;
  end;
end;

% SVD of A = U*S*V'
[U,SA,V] = svd(A);

% Build right hand-side matrix B column-wise
% First build transformation matrix M such that ptilde = p*M
M = zeros(n+1);
for row = 0:n,
  polM = (-S(1,1)-S(1,2)*s)^row*(S(2,1)+S(2,2)*s)^(n-row);
  M(1+row,1:deg(polM)+1) = polM{:};
end;
M = M / sqrt(S(1,2)*S(2,1)-S(1,1)*S(2,2))^n;
% Then build B with components from p'*p-M'*p'*p*M
B = zeros((n+1)*(n+2)/2);
colB = 1;
for row = 1:n+1,
  for col = 1:row,
    mat = zeros(n+1); mat(row,col) = 1; mat(col,row) = 1;
    mat = mat - M'*mat*M;
    B(:, colB) = trivec(mat);
    colB = colB+1;
  end;
end;

% Solve A*h = B*pp = U*SA*V'*h, i.e. h = V*inv(SA)*U'*B*pp
SA = diag(SA(1:n*(n+1)/2,1:n*(n+1)/2)); SA = diag(1./SA); Q = U'*B; 
C = V*SA*Q(1:n*(n+1)/2,:);
if norm(A*C-B,'inf')/norm(B,'inf') > 1e-8,
  error('Inconsistent linear system of equations.');
end;

% end of function BUILDHERMITE
% **************************************************************************

function vec = trivec(mat)
n = size(mat,1);
vec = zeros(n*(n+1)/2,1);
row = 1;
for rowtri = 1:n,
  vec(row:row+rowtri-1) = mat(rowtri, 1:rowtri)';
  row = row+rowtri;
end;

