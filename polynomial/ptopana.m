function stable = ptopana(A,S)
%PTOPANA Robust stability analysis of a polytope of polynomial matrices
%
%   Given a cell array of polynomial matrices A = {A1,A2..,AN},
%   if the instruction
%
%     PTOPANA(A)
%
%   returns 1 then the polytope of polynomial matrices with vertices Ai is
%   robustly stable (i.e. all convex combinations of the Ais are stable).
%   If the instruction returns 0 then we cannot conclude about robust
%   stability.
%
%   With a 2x2 Hermitian matrix S as an additional input argument
%
%     PTOPANA(A,S)
%
%   robust stability is checked within a region of the complex plane
%   {s : S(1,1)+S(1,2)*s+S(2,1)*s'+S(2,2)*s*s' < 0}.

% Robust stability is checked by solving an LMI optimization problem
% following the sufficient condition described in [D. Henrion,
% D. Arzelier, D. Peaucelle, M. Sebek. An LMI condition for robust
% stability of polynomial matrix polytopes. Automatica, Vol. 37, No. 3,
% pp. 461-468, March 2001]. Tolerance for enforcing strict positivity of
% the LMIs is set by parameter PGLOBAL.ZEROING.
%
% The LMI optimization problem is solved with SeDuMi and its
% LMI interface. Both packages must be properly installed.
    
% Author: D. Henrion, January 28, 2002
% Last modified by D. Henrion, May 24, 2002
% Copyright 2002 by PolyX, Ltd.

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

% Polynomial vertices
if isa(A,'double') | isa(A,'pol'), A = {A}; end;
d = deg(A{1}); % degree
n = size(A{1},1); % size
N = length(A); % number of vertices
symb = symbol(A{1}); % variable
for i = 1:N,
 dA = deg(A{i});
 if dA > d, d = dA; end;
 [nA,m] = size(A{1});
 if (nA ~= n) | (nA ~= m),
  error('Inconsistent sizes in input polynomial matrix polytope');
 end;
 newsymb = symbol(A{i});
 if isempty(symb),
   symb = newsymb;
 end;
 if ~isempty(newsymb) & ~strcmp(symb, newsymb),
   error('Inconsistent variable symbols in input polynomial matrix polytope');
 end;
end;

% Stability region
if nargin < 2,
 switch symb,
  case {'z^-1','d'}
   S = [1 0;0 -1];
  case {'z','q'}
   S = [-1 0;0 1];
  otherwise
   S = [0 1;1 0]; 
 end;
end;

% Check stability region
e = eig(S);
if sum(e<=0) == 2,
  error('Stability region is empty.');
elseif sum(e>=0) == 2,
  error('Stability region is the whole complex plane.');
end;

R = [eye(n*d) zeros(n*d,n); zeros(n*d,n) eye(n*d)];
coefA = cell(1,N);
for i = 1:N,
 coefA{i} = zeros(n,(d+1)*n);
 for j = 0:d,
  coefA{i}(:,1+j*n:(j+1)*n) = A{i}{j};
 end;
end;

% LMI variables
lmi = sdmpb('Robust stability analysis');
varP = cell(1,N);
for i = 1:N,
 [lmi, varP{i}] = sdmvar(lmi, d*n, 's', ['P' int2str(i)]);
end;
[lmi, varQ] = sdmvar(lmi, 2*d*n, n, 'Q');

% define LMI
% Pi > 0
for i = 1:N,
 [lmi, index] = sdmlmi(lmi, d*n);
 lmi = sdmineq(lmi, -index, varP{i});
end;

% R'*kron(S,Pi)*R+Ai'*Q'*R+R'*Q*Ai < 0
for i = 1:N,
 [lmi, index] = sdmlmi(lmi, (d+1)*n);
% lmi = sdmineq(lmi, index, varP{i}, [], R, S); % R'*kron(S,Pi)*R
 lmi = sdmineq(lmi, index, varP{i}, R', [], S); % R'*kron(S,Pi)*R
 lmi = sdmineq(lmi, index, varQ, R', coefA{i}); % R'*Q*Ai+()'
end;

% strict LMI
lmi = lmi - tol;

% solve LMI
pars.fid = verbose;
lmi = sdmsol(lmi,pars);

feas = get(lmi, 'feas');
stable = (feas > 0);

% check LMI
%P = cell(1,N);
%M = cell(1,N);
%Q = lmi(varQ);
%for i = 1:N,
% P{i} = lmi(varP{i});
% M{i} = R'*kron(S,P{i})*R+coefA{i}'*Q'*R+R'*Q*coefA{i};
%end;
