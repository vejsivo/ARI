function [F0,F1] = ellides2(An,Ap,B,C,D,S)
% ELLIDES2 - Robust stabilization by PI feedback of an ellipsoid of
%           second-order systems
%
% Given quadratic polynomial matrices A0(s) and A1(s), constant input
% matrix B, constant output matrix C and quadratic polynomial
% matrix D(s), the instruction
%
%   [F0,F1] = ELLIDES2(A0,A1,B,C,D)
%
% attempts to find two constant feedback matrices F0 and F1 such
% that the second order system (A0(s) + Q'*A1(s))*x = B*u, y = C*x
% affected by ellipsoidal uncertainty (uncertain parameter Q is
% such that NORM(Q) <= 1) is robustly stabilized by the
% proportional and derivative feedback u = (F0+F1*s)*y
%
% Polynomial matrix D(s) corresponds to a central (or nominal) closed-loop
% polynomial matrix, for example D(s) = (A0(s)+B*F0*C+B*F1*C*s) 
% obtained by stabilizing the nominal plant.
%
% By default, A1 is zero, B and C are identity
%
% With the optional input argument
%
%   [F0,F1] = ELLIDES2(A0,A1,B,C,D,S)
% 
% one can specify an Hermitian 2x2 stability region matrix S such that
% each closed-loop pole Z must satisfy the quadratic inequality
% S(1,1)+S(1,2)*Z+S(2,1)*Z'+S(2,2)*Z*Z' < 0. Standard choices are
% S=[0 1;1 0] (left half-plane) and S=[-1 0;0 1] (unit disk).
% The default value of S is determined by the variable symbol of input
% polynomials A0,A1 and D.

% Written by Didier Henrion, May 15, 2002.
% Last modified by Didier Henrion, May 24, 2002.
% Copyright 2002 by PolyX, Ltd.  

global PGLOBAL;
  
eval('PGLOBAL.FORMAT;',...
     'error(''Use PINIT to initialize the Polynomial Toolbox.'');');

if nargin < 5,
  error('Invalid number of input arguments');
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

D = pol(D);
n = size(D,1);
if deg(D) > 2, error('Invalid degree of fifth input polynomial'); end;
if size(D,2) ~= n, error('Invalid size of fifth input polynomial'); end;

symb = symbol(D);
An = pol(An); Ap = pol(Ap);
M = {An,Ap};
for i = 1:2,
 newsymb = symbol(M{i});
 if isempty(symb),
  symb = newsymb;
 elseif ~isempty(newsymb) & ~strcmp(symb, newsymb),
  error('Inconsistent variable symbols in input arguments');
 end;
end;

if deg(An) > 2, error('Invalid degree of first input polynomial'); end;
if size(An) ~= [n n], error('Invalid size of first input polynomial'); end;

if deg(Ap) > 2, error('Invalid degree of second input polynomial'); end;
if isempty(Ap), Ap = pol(zeros(n)); end;
if size(Ap,2) ~= n, error('Invalid size of second input polynomial'); end;
q = size(Ap,1); % perturbation size

if isempty(B), B = eye(n); end;
if size(B,1) ~= n, error('Invalid size of third input matrix'); end;
m = size(B,2);

if isempty(C), C = eye(n); end;
if size(C,2) ~= n, error('Invalid size of third input matrix'); end;
p = size(C,1);

if nargin < 6, S = []; end; 
if isempty(S),
 switch symb,
  case {'z^-1','d'}
   S = [1 0;0 -1];
  case {'z','q'}
   S = [-1 0;0 1];
  otherwise
   S = [0 1;1 0]; 
 end;
end;
if ~isa(S,'cell'),
 S = {S};
end;
nS = length(S);
for i = 1:nS,
 if ~isa(S{i},'double') | (size(S{i}) ~= [2 2]),
  error('Invalid stability region');
 end;
end;

lmi = sdmpb;
[lmi, varF0] = sdmvar(lmi, m, p, 'F0');
[lmi, varF1] = sdmvar(lmi, m, p, 'F1');
for j = 1:nS,
 [lmi, varP{j}] = sdmvar(lmi, 2*n, 's');
end;
%[lmi, varrho] = sdmvar(lmi, 1, 1);
[lmi, vargamma] = sdmvar(lmi, 1, 1);
[lmi, varnorm] = sdmvar(lmi, 1, 1);

coefD = zeros(n,3*n);
coefN = zeros(n,3*n);
coefAp = zeros(q,3*n);
for d = 0:2,
 coefD(:,1+d*n:(d+1)*n) = D{d};
 coefN(:,1+d*n:(d+1)*n) = An{d};
 coefAp(:,1+d*n:(d+1)*n) = Ap{d};
end;
for j = 1:nS,
 name = ['Region ' int2str(j)];
 [lmi, index] = sdmlmi(lmi, 3*n+q, name);
 left = [eye(3*n); zeros(q,3*n)];
 right = [eye(3*n) zeros(3*n,q)];
 lmi = sdmineq(lmi, -index, 0, left*coefD', coefN*right);
 lmi = sdmineq(lmi, -index, varF1, left*coefD'*B, ...
	       [zeros(p,n) C zeros(p,n)]*right);
 lmi = sdmineq(lmi, -index, varF0, left*coefD'*B, ...
	       [C zeros(p,2*n)]*right);
 Q = kron([1 0 0;0 1 0;0 1 0;0 0 1],eye(n));
 lmi = sdmineq(lmi, index, varP{j}, left*Q', [], S{j});
 lmi = sdmineq(lmi, index, vargamma, left*coefD');
 left = [zeros(3*n,q); eye(q)];
 right = [eye(3*n) zeros(3*n,q)];
 % lmi = sdmineq(lmi, -index, varrho, left*coefAp, right);
 lmi = sdmineq(lmi, -index, 0, left*coefAp, right);
 lmi = sdmineq(lmi, -index, vargamma, left);
end;

% minimization of norm of feedback matrices
[lmi, index] = sdmlmi(lmi, p+2*m, 'Min norm');
left = [eye(p); zeros(2*m,p)];
lmi = sdmineq(lmi, -index, varnorm, left);
left = [zeros(p,m); eye(m); zeros(m)];
right = [eye(p) zeros(p,2*m)];
lmi = sdmineq(lmi, -index, varF0, left, right);
left = [zeros(p+m,m); eye(m)];
lmi = sdmineq(lmi, -index, varF1, left, right);
term = [zeros(p) zeros(p,2*m); zeros(2*m,p) eye(2*m)];
lmi = sdmineq(lmi, -index, 0, term, term/2);
if verbose,
  disp('Minimize norm of feedback matrices');
end;

lmi = sdmobj(lmi, varnorm, -1, 1);

% maximization of uncertainty level
% lmi = sdmobj(lmi, varrho, 1, 1);

pars.fid = verbose;
lmi = sdmsol(lmi,pars,1e8);

feas = get(lmi, 'feas');

% Retrieve controller

if (feas == 0) | (feas == 1),

  P = cell(nS);
  for j = 1:nS,
    P{j} = lmi(varP{j});
  end;
  F0 = lmi(varF0);
  F1 = lmi(varF1);
  % gamma = lmi(vargamma);
  % rho = lmi(varrho);

 if verbose,
   disp('A robustly stabilizing controller was found');
 end;

else

 F0 = []; F1 = [];
 if verbose,
   disp('No controller was found');
 end;
 
end;
