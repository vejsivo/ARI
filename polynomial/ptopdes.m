function [x,y] = ptopdes(a,b,d,S,pid,gamma)
% PTOPDES - Robust stabilization of a polytope of polynomials
%  
% Given two series of polynomials A = {A1,A2,..,AN} and B = {B1,B2,..,BN}
% and a stable polynomial C, the instruction
%
%  [X,Y] = PTOPDES(A,B,C)
%
% attempts to find two polynomials X and Y such that the proper
% controller Y/X of order M = DEG(C)-DEG(A) robustly stabilizes the whole
% polytope of scalar plants with vertices B{i}/A{i}.
%
% Roughly speaking, polynomial C corresponds to the central (or nominal)
% closed-loop polynomial, for example obtained by stabilizing one vertex
% (see macro AXBYC). It fixes the order of the sought controller. If no
% controller is returned, it does not mean that the polytopic plant is
% not robustly stabilizable.
%
% With the optional input arguments
%
%  [X,Y] = PTOPDES(A,B,C,S,PID,GAMMA)
%
% one can specify an Hermitian 2x2 stability region matrix S such that
% each closed-loop pole Z must satisfy the quadratic inequality
% S(1,1)+S(1,2)*Z+S(2,1)*Z'+S(2,2)*Z*Z' < 0. Standard choices are
% S=[0 1;1 0] (left half-plane) and S=[-1 0;0 1] (unit disk).
% The default value of S is determined by the variable symbol of input
% polynomials A,B,C.
%
% When the fifth input argument PID is present, then the controller
% denominator is set to X = s. So if Y/X has order one, it is
% proportional-integral Y/X = Kp+Ki/s, and if Y/X has order two,
% it is a PID Y/X = Kp+Ki/s+Kd*s.
%  
% The sixth input argument GAMMA is a positive real number used to
% ensure stability via strict positive realness (default value = 
% global zeroing tolerance, see GPROP)

% The function is based on the theory described in [D. Henrion,
% D. Arzelier, D. Peaucelle. Positive Polynomial Matrices and
% Improved LMI Robustness Conditions. Proceedings of the IFAC World
% Congress}, Barcelona, Spain, July 2002.]
%  
% The LMI optimization problem is solved with SeDuMi and its
% LMI interface. Both packages must be properly installed.
 
% Written by Didier Henrion, April 4, 2001.
% Last modified by Didier Henrion, December 5, 2002.
% Copyright 2001-2002 by PolyX, Ltd.  

global PGLOBAL;

eval('PGLOBAL.FORMAT;',...
     'error(''Use PINIT to initialize the Polynomial Toolbox.'');');

if nargin < 3,
  error('Invalid number of input arguments.');
end;

verbose = strcmp(PGLOBAL.VERBOSE, 'yes');
tol = PGLOBAL.ZEROING;

% Check whether SeDuMi and its LMI interface are installed
if exist('sedumi') ~= 2,
  error('SeDuMi is not properly installed.');
end;
if exist('@sdmpb/sdmpb') ~= 2,
  error('LMI interface is not properly installed.');
end;

% Parse input arguments

% Polynomial vertices (1st and 2nd input args)

if isa(a,'double') | isa(a,'pol'), a = {a}; end;
if isa(b,'double') | isa(b,'pol'), b = {b}; end;

N = length(a); % number of vertices
n = -Inf; % system order
symb = []; % variable symbol
for i = 1:N,
 a{i} = pol(a{i}); b{i} = pol(b{i});
 da = deg(a{i}); db = deg(b{i});
 newsymb = symbol(a{i});
 if isempty(symb),
   symb = newsymb;
 end;
 if ~isempty(newsymb) & ~strcmp(symb, newsymb),
   error('Inconsistent variable symbols in input arguments.');
 end;
 newsymb = symbol(b{i});
 if isempty(symb),
   symb = newsymb;
 end;
 if ~isempty(newsymb) & ~strcmp(symb, newsymb),
   error('Inconsistent variable symbols in input arguments.');
 end;
 newsymb = symbol(d);
 if isempty(symb),
   symb = newsymb;
 end;
 if ~isempty(newsymb) & ~strcmp(symb, newsymb),
   error('Inconsistent variable symbols in input arguments.');
 end;
 if da > n, n = da; end; % system order
 if db > n, n = db; end; % system order 
 % make denominator monic
 if abs(lcoef(a{i})) > 0,
  if abs(lcoef(a{i})) < tol,
   warning(['Small leading coefficient in denominator polynomial. Results' ...
	  ' may be inaccurate']);
  end;
  b{i} = pol(b{i} / lcoef(a{i})); 
  a{i} = pol(a{i} / lcoef(a{i}));
 end;
end;

% Central polynomial (3rd input arg)

% make characteristic polynomial monic
if abs(lcoef(d)) > 0,
 if abs(lcoef(d)) < tol,
  warning(['Small leading coefficient in central polynomial. Results' ...
	  ' may be inaccurate']);
 end;
 d = pol(d / lcoef(d));
end;
% controller order
m = deg(d)-n;
if m < 0, error(['Degree of central polynomial less than degree of' ...
		 ' plant.']); end;
p = n+m; % order of SPR function

% Stability region (4th input arg)

if nargin < 4,
 S = [];
end;
if ~isempty(S),
 if ~isa(S,'double'),
   error('Invalid fourth input argument.');
 end;
 S = (S+S')/2;
else,
  % stability region
 switch symb,
  case {'s','p'}
   S = [0 1;1 0]; 
  case {'z^-1','d'}
   S = [1 0;0 -1];
  case {'z','q'}
   S = [-1 0;0 1];
  otherwise
   error('Invalid input polynomials.');
 end;
end;

% Structure constraints (5th input arg)

if nargin > 4,
 struct = ~isempty(pid);
else
 struct = false;
end;
if struct
 struct = [zeros(1,m+1) Inf*ones(1,m+1)];
 if m == 0,
  struct(1) = 1;
 else
  struct(2) = 1;
 end;
else
 struct = [Inf*ones(1,m) 1 Inf*ones(1,m+1)];
end;

% Margin for ensuring SPRness

if nargin == 6,
  if ~isa(gamma,'double') | (gamma < 0),
    error('Invalid sixth input argument: must be a positive real');
  end;
else
  gamma = tol;
end;

% Check stability of characteristic polynomial

rd = roots(d);
for i = 1:length(rd),
  if S(1,1)+S(2,1)*rd(i)+S(1,2)*rd(i)'+S(2,2)*rd(i)*rd(i)' >= 0,
    error('Unstable central polynomial.');
  end;
end;
    
if verbose,
 disp('Build LMI.');
end;

% Build projection matrices H{k} to check SPRness condition
% along the stability boundary

H = projection(S,p);

% Build vector coefficients, vertex Sylvester matrices

dcoef = d{:}';
sylv = cell(N,1);
for i = 1:N,
  acoef = a{i}{:}';
  bcoef = b{i}{:}';
  sylv{i} = zeros(n+m+1,2*(m+1));
  for j = 0:m,
    sylv{i}(1+j:length(acoef)+j,1+j) = acoef;
    sylv{i}(1+j:length(bcoef)+j,2+m+j) = bcoef;
  end;
end;

% Given d(s), seek the polynomials x(s),y(s) such that each vertex
% polynomial ci(s)=ai(s)*x(s)+bi(s)*y(s) makes the polynomial
% ci'(s)*d(s)+ci(s)*d'(s)-gamma*d'(s)*d(s) strictly positive real
% along the stability boundary a+b*(s+s')+c*s*s'=0 for a given positive
% (but arbitrarily small) scalar gamma

% Build LMI with decision variables:
% gamma (1) positive;
% r (1), x (m+1), y (m+1) in a cone s.t. norm([x;y])<=r
% X1 (p+1)^2, .., XN (p+1)^2 SDP

K.l = 1;
K.q = 1+2*(m+1);
K.s = (p+1)*ones(1,N);
nvar = K.l+K.q+sum(K.s.^2);
nsc = sum(~isinf(struct)); % structure constraints

% p+1 constraints trace(Hk*Xi)=trace(Hk*(ci*d'+d*ci'-gamma*d*d'))
% for each vertex i, where ci = sylvester(ai,bi)*[x;y]

A = sparse(N*(p+1)+1+nsc,nvar);
B = sparse(size(A,1),1); C = sparse(size(A,2),1);
for i = 1:N, % for each vertex
  for k = 1:p+1, % for each component in the positive polynomial
    row = (i-1)*(p+1)+k;
    A(row,1) = -dcoef'*H{k}*dcoef; % times gamma
    A(row,3:2+2*(m+1)) = 2*dcoef'*H{k}*sylv{i}; % times [x;y]
    col = 3+2*(m+1)+(i-1)*(p+1)^2;
    A(row,col:col+(p+1)^2-1) = -H{k}(:)';
  end;
end;

A(N*(p+1)+1,1) = 1; B(N*(p+1)+1) = gamma; % gamma

% additional structural constraints
row = N*(p+1)+2;
for i = 1:2*(m+1),
  if ~isinf(struct(i)),
   A(row,2+i) = 1; B(row) = struct(i); row = row + 1;
  end;
end;

C(2) = 1; % minimize Euclidean norm of controller vector [x;y]

% Solve LMI

if verbose,
 disp('Solve LMI.');
end;

if verbose, pars.fid = 1; else pars.fid = 0; end;

[X,Y,info] = sedumi(A,B,C,K,pars);

feas = 1;
if any([info.pinf info.dinf info.numerr]),
  if abs(1-info.feasratio) < 1e-1,
    if verbose,
      disp('Problem is almost feasible.');
    end;
  else  
    if verbose,
     disp('No robustly stabilizing controller was found.');
    end;
    x = []; y = []; feas = 0;
  end;
end;

% Retrieve controller

if feas,
  
  gamma = X(1); r = X(2);
  if verbose,
   disp(['SPRness gamma = ' num2str(gamma) ' (must be >0).']);
   disp(['Controller norm = ' num2str(r) ' (must be small).']);
  end;
  
  for i = 1:N,
    col = 3+2*(m+1)+(i-1)*(p+1)^2;
    P = zeros(p+1);
    P(:) = X(col:col+(p+1)^2-1);
    mineigX = min(eig(P));
    if mineigX < 0,
     error(['Min eig X(' int2str(i) ') = ' num2str(mineigX) ...
	  ' (must be >0).']);
    end;
  end
  
  % Retrieve controller polynomials

  x = pol(X(3:m+3)',m,symb);
  y = pol(X(m+4:2*m+4)',m,symb);
  
  % Enforce structure constraints
  
  for i = 0:m,
    if ~isinf(struct(i+1)), x{i} = struct(i+1); end;
    if ~isinf(struct(i+m+2)), y{i} = struct(i+m+2); end;
  end;
    
  problem = 0;
  for i = 1:N,
    r = roots(a{i}*x+b{i}*y);
    maxr = -Inf;
    for j = 1:length(r),
      maxr = max(maxr, S(1,1)+S(2,1)*(r(j)+r(j)')+S(2,2)*r(j)*r(j)');
    end;
    if maxr > 0,
     warning(['Max root(' int2str(i) ') = ' num2str(maxr) ...
	    ' (must be <0).']);
     disp(['Try to move the roots of central polynomial away from the' ...
	    ' stability boundary or try to increase parameter gamma']);
     problem = 1;
    end;
  end;
  if problem
    x = []; y = [];
  end;
  
end;

function H = projection(S,n)
% Given a stability region matrix S = [a b;b' c] and a degree n,
% the command
%
%  H = PROJECTION(S,n)
%
% builds projection matrices H{1},..,H{n+1} of size n+1 to ensure
% positivity along the curve a+b*s+b'*s'+c*s*s'=0 of a polynomial
% p(s) = [1 s s^2 ..s^n]*p. More specifically, p(s) >= 0 along the curve
% if and only if there exists a positive semidefinite matrix X of size
% n+1 satisfying the n+1 scalar constraints p'*H{i}*p = TRACE(H{i}*X).

% vectorize basis matrices
B = zeros(n*(n+1)/2,(n+1)*(n+2)/2);
rowB = 1;
for row = 1:n,
  for col = 1:row,
   M = zeros(n+1);
   M(row:row+1,col:col+1) = S;
   if row ~= col,
     M(col:col+1,row:row+1) = M(col:col+1,row:row+1)+S;
   end;
   B(rowB,:) = trivec(M)'; rowB = rowB+1;
  end;
end;
% extract null space
N = null(B);
if size(N,2) ~= n+1,
  error('Invalid null space for projection matrix.');
end;
% N = rref(N')'; % reduced echelon form (for nice display)

% retrieve projection matrices
H = cell(n+1,1);
for col = 1:n+1,
  H{col} = trimat(N(:,col));
end;

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
  mat(rowtri, 1:rowtri) = vec(row:row+rowtri-1)';
  row = row+rowtri;
end;
mat = (mat+mat')/2;

