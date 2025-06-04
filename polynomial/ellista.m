function [P,pc] = ellista(pc,S,center)
%ELLISTA Ellipsoidal approximation of the stability domain
%
%   Given a stable monic polynomial P0 of degree N, the instruction
%
%    [P,PC] = ELLISTA(P0)
%
%   returns a positive definite matrix P and a polynomial PC such that the
%   ellipsoid E = {p: (p-pc)*P*(p-pc) <= 1} is an inner approximation of
%   the stability domain in the space of coefficients of monic
%   polynomials of degree N, i.e. such that every vector
%   p = [p0;p1;p2..;1] in E parametrizes a stable polynomial
%   P = p0+p1*s+..+s^N. Ellipsoid E contains the vector corresponding to
%   the given polynomial P0. The trace of P is minimized so that the
%   volume of E is indirectly maximized.
%
%   With a 2x2 Hermitian matrix S as an additional input argument,
%   the instruction
%
%     [P,PC] = ELLISTA(P0,S)
%
%   returns an ellipsoid E such that every vector in E parametrizes a
%   polynomial with its roots in the region of the complex plane
%   S = {s : S(1,1)+S(1,2)*s+S(2,1)*s'+S(2,2)*s*s' < 0}. If matrix S is not
%   specified, then the default stability region corresponds
%   to the variable symbol in input polynomial P0.
%
%   With the syntax
%
%     P = ELLISTA(PC,S,'center')
%
%   the monic polynomial center PC of ellipsoid E is specified
%   instead of P0.
  
% The ellipsoid is obtained by solving an LMI optimization problem
% built after applying the S-procedure to relax the non-convex BMI
% corresponding to the Hermite stability criterion, as explained in
% [D. Henrion, D. Peaucelle, D. Arzelier, M. Sebek. Ellipsoidal
% Approximation of the Stability Domain of a Polynomial. Proceedings of
% the European Control Conference, pp. 384-389, Porto, Portugal, 2001]
% Note that the LMI formulation slightly differs from the one found
% in the above reference. In the current implementation it is less
% conservative.
%
% The LMI optimization problem is solved with SeDuMi and its
% LMI interface. Both packages must be properly installed.
 
% Author: D. Henrion, September 17, 2001.
% Last modified by D. Henrion, May 24, 2002.
% Copyright 2001-2002 by PolyX, Ltd.

global PGLOBAL;

eval('PGLOBAL.FORMAT;',...
     'error(''Use PINIT to initialize the Polynomial Toolbox.'');');

verbose = strcmp(PGLOBAL.VERBOSE, 'yes');
tol = PGLOBAL.ZEROING;

% Check whether SeDuMi and its LMI interface are installed
if exist('sedumi') ~= 2,
  error('SeDuMi is not properly installed.');
end;
if exist('@sdmpb/sdmpb') ~= 2,
  error('LMI interface is not properly installed.');
end;
if exist('@sdmpb/sdmlme') ~= 2,
  error('LMI interface version 1.03 or higher must be installed.');
end;

if nargin < 2,
 S = [];
end;

if isempty(S),
 switch symbol(pc),
  case {'z^-1','d'}
   S = [1 0;0 -1];
  case {'z','q'}
   S = [-1 0;0 1];
  otherwise
   S = [0 1;1 0]; 
 end;
end;

if nargin < 3,
  center = 0;
  if verbose,
    disp('Compute stability ellipsoid with given point.');
  end;
else
  center = 1;
  if verbose,
    disp('Compute stability ellipsoid with given center.');
  end;
end;

% Check stability region
e = eig(S);
if sum(e<=0) == 2,
  error('Stability region is empty.');
elseif sum(e>=0) == 2,
  error('Stability region is the whole complex plane.');
end;

% Is pc monic ?
if abs(lcoef(pc)-1) > tol,
  error('Input polynomial must be monic.');
end;

% Check roots of pc
pc = pol(pc);
r = roots(pc);
for i = 1:length(r),
  rs = r(i);
  if S(1,1)+S(1,2)*rs+S(2,1)*rs'+S(2,2)*rs*rs' >= 0,
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
 disp('Build Hermite matrix.');
end;
C = hermcoef(S,n);

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

% Build the LMI kron(D,eye(n+1))*H >= kron(eye(n),PP)+G
% with kron(D,eye(n+1))*H = H*kron(D,eye(n+1))
% and D = D' > 0 and G = [Gij], Gij = -Gij'

if verbose,
  disp('Build LMI.');
end;
lmi = sdmpb('Ellipsoidal approximation');

% Decision variables: D, P (or PP) and G
[lmi, Dindex] = sdmvar(lmi, n, 's', 'D');
if center,
 [lmi, Pindex] = sdmvar(lmi, n, 's', 'P');
else
 [lmi, PPindex] = sdmvar(lmi, n+1, 's', 'PP');
end;

Gindex = zeros(n*(n-1)/2,1);
for i = 1:n*(n-1)/2,
  [lmi, index] = sdmvar(lmi, n+1, 'as', ['G' int2str(i)]); Gindex(i) = index;
end;
Gvar = zeros((n+1)*n);
ind = 1;
for row = 2:n,
  for col = 1:row-1,
    var = lmi{Gindex(ind)};
    Gvar(1+(row-1)*(n+1):row*(n+1),1+(col-1)*(n+1):col*(n+1)) = var;
    Gvar(1+(col-1)*(n+1):col*(n+1),1+(row-1)*(n+1):row*(n+1)) = var';
    ind = ind + 1;
  end;
end;
[lmi, Gindex] = sdmvar(lmi, Gvar, 'st', 'G');

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

% D > 0
[lmi, index] = sdmlmi(lmi, n, 'D > 0');
lmi = sdmineq(lmi, -index, Dindex);

% kron(D,eye(n+1))*H >= kron(eye(n),PP)+G
[lmi, index] = sdmlmi(lmi, n*(n+1), ...
		      'kron(D,eye(n+1)))*H >= kron(eye(n),PP)+G');
lmi = sdmineq(lmi, index, Gindex);
lmi = sdmineq(lmi, -index, Dindex, 0.5, H, eye(n+1), -1);

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

% commutation kron(D,eye(n+1))*H = H*kron(D,eye(n+1))
[lmi, index] = sdmlme(lmi, n*(n+1), n*(n+1), 'commutation');
lmi = sdmeq(lmi, -index, Dindex, 1, H, eye(n+1), -1);
lmi = sdmeq(lmi, index, Dindex, H, 1, eye(n+1), -1);

% Solve the LMI
if verbose,
  disp('Solve LMI.');
end;
pars.fid = verbose;
lmi = sdmsol(lmi,pars);

feas = get(lmi, 'feas');
if feas >= 0,
 if center,
  D = lmi(Dindex); P = lmi(Pindex); G = lmi(Gindex);
  PP = [-P P*pcoef;pcoef'*P 1-pcoef'*P*pcoef];
  M = kron(D,eye(n+1))*H - kron(eye(n),PP)-G;
 else
  D = lmi(Dindex); PP = lmi(PPindex); G = lmi(Gindex);
  M = kron(D,eye(n+1))*H - kron(eye(n),PP) -G;
  P11 = PP(1:n,1:n); P12 = PP(1:n,n+1); P22 = PP(n+1,n+1);
  pcoef = -P11\P12;
  P = -P11/(P22+P12'*pcoef);
  pc = pol([pcoef' 1],n,symbol(pc));
 end;
 if verbose,
   disp(['Norm D=' num2str(norm(D)) ', norm G=' ...
	 num2str(norm(G)) ', logdet P=' num2str(log(det(P)))]);
 end;
else
 P = []; pc = [];
end;

% end of function ELLISTA
% **************************************************************************

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

