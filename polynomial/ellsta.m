function [P,pc] = ellsta(a,b,c,pc,center)
%ELLSTA Ellipsoidal approximation of the stability domain
%
%   Given a stability region S = {s : a+b*s+b'*s'+c*s*s' < 0} in the
%   complex plane and a polynomial p0 of degree N with all its roots in S,
%   the instruction
%
%     [P,PC] = ELLSTA(a,b,c,P0)
%
%   returns a positive definite matrix P and a polynomial PC such that the
%   ellipsoid E = {p: (p-pc)*P*(p-pc) <= 1} is an inner approximation of
%   the stability domain in the space of coefficients of monic
%   polynomials of degree N, i.e. such that every vector
%   p = [p0;p1;p2..;1] in E parametrizes a polynomial P = p0+p1*s+..+s^N
%   with all its roots in S. Ellipsoid E contains the vector
%   corresponding to the given polynomial P0. The trace of P is minimized
%   so that the volume of E is indirectly maximized.
%
%   With the syntax
%
%     P = ELLSTA(a,b,c,PC,'center')
%
%   the center PC of ellipsoid E is specified.
  
% The ellipsoid is obtained by solving an LMI optimization problem
% built after applying the S-procedure to relax the non-convex BMI
% corresponding to the Hermite stability criterion, as explained in
% [D. Henrion, D. Peaucelle, D. Arzelier, M. Sebek. Ellipsoidal
% Approximation of the Stability Domain of a Polynomial. Proceedings of
% the European Control Conference, pp. 384-389, Porto, Portugal, 2001]
 
% Author: D. Henrion, September 17, 2001.
% Copyright 1998-2001 by PolyX, Ltd.

global PGLOBAL;

verbose = strcmp(PGLOBAL.VERBOSE, 'yes');

if nargin < 4,
  error('Invalid number of input arguments.');
end;

if nargin < 5,
  center = 0;
  if verbose,
    disp('ELLSTA: Compute stability ellipsoid with given point.');
  end;
else
  center = 1;
  if verbose,
    disp('ELLSTA: Compute stability ellipsoid with given center.');
  end;
end;

% Check stability region
e = eig([a b;b' c]);
if sum(e<=0) == 2,
  error('Stability region is empty.');
elseif sum(e>=0) == 2,
  error('Stability region is the whole complex plane.');
end;

% Check roots of pc
pc = pol(pc);
r = roots(pc);
for i = 1:length(r),
  rs = r(i);
  if a+b*rs+(b*rs)'+c*rs*rs' >= 0,
    error('Roots of input polynomial not in stability region.');
  end;
end;

% Make pc monic
n = deg(pc);
if n < 2,
  error('Degree of input polynomial is less than two.');
end;
pcoef = zeros(n,1);
for i = 1:n,
 pcoef(i) = pc{i-1} / pc{n};
end;

% Build coefficients of the Hermite matrix
if verbose,
 disp('ELLSTA: Build Hermite matrix.');
end;
C = buildhermite(a,b,c,n);

% Build the big matrix H
H = zeros(n*(n+1));
rowtri = 1;
for row = 1:n,
  for col = 1:row,
   mat = trimat(C(rowtri,:));
   H(1+(row-1)*(n+1):row*(n+1),1+(col-1)*(n+1):col*(n+1)) = mat;
   H(1+(col-1)*(n+1):col*(n+1),1+(row-1)*(n+1):row*(n+1)) = mat;
   rowtri = rowtri+1;
  end;
end;

% Build the LMI l*H >= kron(eye(n),PP)+S with S = [Sij], Sij = -Sij'

if verbose,
  disp('ELLSTA: Build LMI.');
end;
lmi = sdmpb('Ellipsoidal approximation');

% Decision variables: l, P (or PP) and S
[lmi, lindex] = sdmvar(lmi, 1, 1, 'lambda');
if center,
 [lmi, Pindex] = sdmvar(lmi, n, 's', 'P');
else
 [lmi, PPindex] = sdmvar(lmi, n+1, 's', 'PP');
end;
Sindex = zeros(n*(n-1)/2,1);
for i = 1:n*(n-1)/2,
  [lmi, index] = sdmvar(lmi, n+1, 'as', ['S' int2str(i)]); Sindex(i) = index;
end;
Svar = zeros((n+1)*n);
ind = 1;
for row = 2:n,
  for col = 1:row-1,
    var = lmi{Sindex(ind)};
    Svar(1+(row-1)*(n+1):row*(n+1),1+(col-1)*(n+1):col*(n+1)) = var;
    Svar(1+(col-1)*(n+1):col*(n+1),1+(row-1)*(n+1):row*(n+1)) = var';
    ind = ind + 1;
  end;
end;
[lmi, Sindex] = sdmvar(lmi, Svar, 'st', 'S');

if center,

 % P > 0
 [lmi, index] = sdmlmi(lmi, n, 'P>0');
 lmi = sdmineq(lmi, -index, Pindex);

else

 % P11 < 0
 [lmi, index] = sdmlmi(lmi, n, 'P11<0');
 lmi = sdmineq(lmi, index, PPindex, [eye(n) zeros(n,1)]);

 % pc belongs to the ellipsoid, i.e. [pc;1]'*PP*[pc;1] >= 0
 % we normalize the constraint to [pc;1]'*PP*[pc;1] = 1
 [lmi, index] = sdmlmi(lmi, 1, '[pc;1]''*PP*[pc;1] >= 1');
 lmi = sdmineq(lmi, -index, PPindex, [pcoef;1]');
 lmi = sdmineq(lmi, index, 0, 1);
 [lmi, index] = sdmlmi(lmi, 1, '[pc;1]''*PP*[pc;1] <= 1');
 lmi = sdmineq(lmi, index, PPindex, [pcoef;1]');
 lmi = sdmineq(lmi, -index, 0, 1);
  
end;

% l*H >= kron(eye(n),PP)+S
[lmi, index] = sdmlmi(lmi, n*(n+1), 'l*H >= kron(eye(n),PP)+S');
lmi = sdmineq(lmi, index, Sindex);
lmi = sdmineq(lmi, -index, lindex, H, 0.5);

if center,
  
 % with PP = [-P P*pc;pc'*P 1-pc'*P*pc], we have kron(eye(n),PP) =
 % -kron(eye(n),R')*kron(eye(n),P)*kron(eye(n),R) + kron(eye(n),Q)
 % for R = [eye(n) -pc], Q = [zeros(n) zeros(n,1); zeros(1,n) 1]
 R = [eye(n) -pcoef]; Q = [zeros(n) zeros(n,1); zeros(1,n) 1];
 lmi = sdmineq(lmi, -index, Pindex, kron(eye(n),R'), kron(eye(n),R), ...
	      0.5*eye(n));
 lmi = sdmineq(lmi, index, 0, kron(eye(n),Q), 0.5);

else
  
 lmi = sdmineq(lmi, index, PPindex, [], [], eye(n));

end;

if center,

 % criterion: min trace P
 lmi = sdmobj(lmi, Pindex, 'tr', -1);

else

 % criterion: max trace P11
 for i = 1:n,
  e = zeros(n+1,1); e(i) = 1;
  lmi = sdmobj(lmi, PPindex, e', e);
 end;

end;

% Solve the LMI
if verbose,
  disp('ELLSTA: Solve LMI.');
end;
pars.fid = verbose;
lmi = sdmsol(lmi,pars);

feas = get(lmi, 'feas');
if feas >= 0,
  if center,
  l = lmi(lindex); P = lmi(Pindex); S = lmi(Sindex);
  PP = [-P P*pcoef;pcoef'*P 1-pcoef'*P*pcoef];
  M = l*H - kron(eye(n),PP)-S;
 else
  l = lmi(lindex); PP = lmi(PPindex); S = lmi(Sindex);
  M = l*H - kron(eye(n),PP)-S;
  P11 = PP(1:n,1:n); P12 = PP(1:n,n+1); P22 = PP(n+1,n+1);
  pcoef = -P11\P12;
  P = -P11/(P22+P12'*pcoef);
  pc = pol([pcoef' 1],n,symbol(pc));
 end;
 if verbose,
   disp(['ELLSTA: lambda=' num2str(l) ', norm S=' ...
	 num2str(norm(S)) ', logdet P=' num2str(log(det(P)))]);
 end;
else
 P = []; pc = [];
end;

% end of function ELLSTA
% **************************************************************************

function C = buildhermite(a,b,c,n)
% BUILDHERMITE: build coefficients of Hermite matrix of a polynomial
% a,b,c: stability region {s: a+b*(s+s')+c*s*s' < 0}
% n: degree of polynomial p(s) = p0+p1*s+p2*s^2+..
% C: coefficients of n-by-n Hermite matrix H = [Hij] of p(s) satisfying
%    [H11;H12;H22;H13;H23;H33..] = C*[p0*p0;p1*p0;p1*p1;p2*p0;p2*p1;p2*p2..]
  
S = [a b;b c]; % stability region
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
[U,S,V] = svd(A);

% Build right hand-side matrix B column-wise
% First build transformation matrix M such that ptilde = p*M
M = zeros(n+1);
for row = 0:n,
  polM = (-a-b*s)^row*(b+c*s)^(n-row); M(1+row,1:deg(polM)+1) = polM{:};
end;
M = M / sqrt(b^2-a*c)^n;
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

% Solve A*h = B*pp = U*S*V'*h, i.e. h = V*inv(S)*U'*B*pp
S = diag(S(1:n*(n+1)/2,1:n*(n+1)/2)); S = diag(1./S); Q = U'*B; 
C = V*S*Q(1:n*(n+1)/2,:);
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

function mat = trimat(vec)
m = length(vec);
n = (1+sqrt(1+8*m))/2-1;
mat = zeros(n);
row = 1;
for rowtri = 1:n,
  mat(rowtri, 1:rowtri) = vec(row:row+rowtri-1);
  row = row+rowtri;
end;
mat = (mat+mat')/2;

