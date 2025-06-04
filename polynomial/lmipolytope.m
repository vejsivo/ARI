function s = lmipolytope(A,option)
%LMIPOLYTOPE Robust stability analysis of a polytope of polynomial
%            matrices via LMIs
%
% The command
%
%    S = LMIPOLYTOPE({A1 A2 .. AN})
%
% checks stability of the polytope of polynomial matrices whose vertices
% are A1, A2, .., AN. All the vertices Ai must be square nonsingular
% polynomial matrices of the same dimension and degree and no zero at infinity.
%
% Stability is checked via an LMI feasibility problem
%
% If S = 1 then it means that the polytope of polynomial matrices is stable.
% If S = 0 then it means that the macro failed to conclude about stability
% of the polytope of polynomial matrices. In particular, it does not
% necessarily mean that the polytope is not stable.
%
% By default, stability is checked with a parameter-dependent Lyapunov matrix.
% However, with the syntax
%
%    T = LMIPOLYTOPE({A1 A2 .. AN}, 'quad')
%
% the macro seeks a unique (parameter-independent) Lyapunov matrix over
% the whole polytope, which may be less time-consuming but more conservative.
% In other words, T = 1 implies S = 1 and S = 0 implies T = 0.

% Author: Didier Henrion, January 27, 2000.
% Updated to 3.0 by D. Henrion, September 1, 2000.
% Modified by J. Jezek, August 24, 2001,  arg checking
% Copyright 2000 by Polyx, Ltd.

global PGLOBAL;
eval('PGLOBAL.VERBOSE;', 'painit;');

if nargin<1
   error('Not enough input arguments.');
end;

quadratic = 0;
if nargin == 2,
 quadratic = strcmp(option, 'quad');
end;

if ~isa(A, 'cell'),
 A = {A};
end;

[rA,cA] = size(A);
N = rA*cA; % number of vertices
if N == 0
   error('No vertices.');
end;

% parse input vertices (size, degree, symbol)

for i = 1:N,
 eval('A{i} = pol(A{i});', ...
      'error([''Argument number '' int2str(i) '' is not polynomial matrix.'']);');
 if i == 1,
  [n,m] = size(A{i});
  if (n ~= m),
   error('Vertex number 1 is not square.');
  end;
  if issingular(A{i}),
   error('Vertex number 1 is singular.');
  end;
  d = deg(A{i});
  sym = symbol(A{i});
 else
  [n2,m2] = size(A{i});
  if (n2 ~= m2),
   error(['Vertex number ' int2str(i) ' is not square.']);
  elseif (n2 ~= n),
   error(['Vertex number ' int2str(i) ' should be ' int2str(n) ...
         'x' int2str(n) '.']);
  elseif issingular(A{i}),
   error(['Vertex number ' int2str(i) ' is singular.']);
  end;
  d = max(d, deg(A{i}));
  if isempty(sym),
   sym = symbol(A{i});
  elseif ~strcmp(symbol(A{i}), sym),
   warning(['Vertex number ' int2str(i) ' has inconsistent symbol.']);
   disp(['Symbol is assumed to be [' sym '].']);
   arg = A{i}; symbol(arg,sym); A{i} = arg;
  end;
 end;
end;

verbose = strcmp(PGLOBAL.VERBOSE, 'yes');

if isempty(d) | d == 0, % only  empty vertices  or  constant vertices
 s = logical(1); return;
end;

% stability region

if strcmp(sym,'s') | strcmp(sym,'p'),
   % closed right half plane
   B = [0 1;1 0];
elseif strcmp(sym,'z^-1') | strcmp(sym,'d'),
   % closed unit disk
   B = [1 0;0 -1];
elseif strcmp(sym,'z') | strcmp(sym,'q'),
   % closed outside of unit disk
   B = [-1 0;0 1];
end;

if N == 1,
  if verbose,
    disp('LMIPOLYTOPE: Polytope has only one vertex: call LMIANALYSIS.');
  end;
  s = lmianalysis(A{1});
  return;
end;
  
% retrieve matrix coefficients

AA = cell(N,1);
for i = 1:N,
 V = A{i}; dV = max(deg(V), 0);
 AA{i} = [V{:} zeros(n, (d-dV)*n)];
end;

% build LMI

PP = [zeros(d*n,n) eye(d*n)];
QQ = [eye(d*n) zeros(d*n,n)];

if verbose,
 disp('LMIPOLYTOPE: Build LMI.');
end;

setlmis([]);

if (quadratic == 0),

 if verbose,
  disp('LMIPOLYTOPE: Seek parameter-dependent Lyapunov matrix.');
 end;

 % parameter-dependent Lyapunov matrix

 varF = lmivar(2,[d*n n]);
 varG = lmivar(2,[d*n n]);

 for i = 1:N,
  varP{i} = lmivar(1,[d*n 1]);
  LMI = newlmi;
  lmiterm([-LMI 1 1 varP{i}], 1, 1);
 end;

 for i = 1:N,
  LMI = newlmi;
  lmiterm([LMI 1 1 varP{i}], B(1,1), 1);
  lmiterm([LMI 2 1 varP{i}], B(2,1), 1);
  lmiterm([LMI 2 2 varP{i}], B(2,2), 1);
  lmiterm([LMI 1 3 varF], 1, 1);
  lmiterm([LMI 2 3 varG], 1, 1);
  lmiterm([LMI 0 0 0], [QQ; PP; AA{i}]);
 end;

else

 % quadratic stability

 if verbose,
  disp('LMIPOLYTOPE: Seek unique Lyapunov matrix.');
 end;
 
 varF = lmivar(2,[d*n n]);
 varG = lmivar(2,[d*n n]);

 varP = lmivar(1,[d*n 1]);
 LMI = newlmi;
 lmiterm([-LMI 1 1 varP], 1, 1);

 for i = 1:N,
  LMI = newlmi;
  lmiterm([LMI 1 1 varP], B(1,1), 1);
  lmiterm([LMI 2 1 varP], B(2,1), 1);
  lmiterm([LMI 2 2 varP], B(2,2), 1);
  lmiterm([LMI 1 3 varF], 1, 1);
  lmiterm([LMI 2 3 varG], 1, 1);
  lmiterm([LMI 0 0 0], [QQ; PP; AA{i}]);
 end;

end;

LMISYS = getlmis;

if verbose,
 disp(['LMIPOLYTOPE: Number of decision variables = ' ...
  int2str(decnbr(LMISYS)) '.']);
end;

% solve LMI

if verbose,
 disp('LMIPOLYTOPE: Solve LMI feasibility problem.');
end;

options(5) = 1;
[tmin,xfeas] = feasp(LMISYS, options);

if isempty(tmin),
 error('LMI is infeasible !');
end;

s = (tmin <= 0);

if verbose,
 if s,
  disp('LMIPOLYTOPE: Polytope is stable.');
 else
  disp('LMIPOLYTOPE: Cannot conclude about stability of polytope.');
 end;
end;

%end .. lmipolytope
