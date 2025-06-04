function [nc,dc] = lmisimstab(pn,pd,d,tol)
%LMISIMSTAB Simultaneous stabilization of a family of SISO systems via LMIs
%
% The command
%
%   [NC,DC] = LMISIMSTAB(N,D,ORDER)
%
% attempts to find a continuous-time SISO controller NC/DC of order less than
% or equal to ORDER simultaneously stabilizing a series of SISO linear systems
% N{i}/D{i}, i.e. such that characteristic polynomials NC*N{i}+DC*D{i} are all 
% Hurwitz stable. If the third argument is unspecified, then ORDER is
% assumed to be equal to the order of the systems.
%
% When there are 3 or more plants to be simultaneously stabilized, the
% design problem is approached via a non-convex rank-constrained LMI problem,
% in turned approached with macro LMIRANK. This macro is based on a
% cone-complementarity heuristic without guarantee of convergence. 
%
% Therefore, if NC and DC are empty, then it means that the macro failed to 
% find a simultaneously stabilizing controller. However, it does mean that
% such a controller does not exist.
%
% A relative accuracy TOL can be specified as a fourth input argument.
% It is used in macro LMIRANK. Its default value is the global zeroing
% tolerance.

% The rank-constrained LMI formulation of the simultaneous stabilization
% problem is described in D. Henrion, S. Tarbouriech, M. Sebek "Rank-one
% LMI Approach to Simultaneous Stabilization of Linear Systems" Systems and
% Control Letters, Vol. 38, pp. 79-89, 1999.

% Author: Didier Henrion, February 8, 2000.
% Modified by Didier Henrion, March 20, 2000.
% Updated to 3.0 by Didier Henrion, September 1, 2000.
% Modified by Jan Jezek, Aug 2001, arg checking
% Copyright 2000 by Polyx, Ltd.

global PGLOBAL;
eval('PGLOBAL.ZEROING;', 'painit;');

if nargin < 2,
 error('Not enough input arguments.');
elseif nargin < 3,
 d = [];
 tol = [];
elseif nargin < 4,
 tol = [];
end;

if isempty(tol),
   tol = PGLOBAL.ZEROING;
elseif ~isa(tol,'double') | length(tol)~=1 | ...
      ~isreal(tol) | tol<0 | tol>1,
   error('Invalid tolerance.');
end;

verbose = strcmp(PGLOBAL.VERBOSE, 'yes');

if ~isa(pn, 'cell'),
 pn = {pn}; % numerator polynomials
end;
if ~isa(pd, 'cell'),
 pd = {pd}; % denominator polynomials
end;

N = length(pn); % number of systems to be stabilized
Nn = length(pd);

if N ~= Nn,
 error('1st and 2nd input arguments have inconsistent dimensions.');
end;

if N <= 1,
 error('At least two systems must be specified.');
end;

% parse SISO systems (size, degree, symbol)

inconsym = 0;
for i = 1:N,

   % numerators
   eval('arg = pol(pn{i});', 'error(peel(lasterr));');
   if size(arg) ~= [1 1],
      error('Numerator polynomials must be scalars.');
   end;
   if ~(isempty(symbol(arg)) | strcmp(symbol(arg),'s') | ...
         strcmp(symbol(arg),'p')),
      error('Continuous time only.');
   end;
   if i == 1,
      ds = deg(arg);
      sym = symbol(arg);
      pn{i} = arg;
   else
      ds = max(ds, deg(arg));
      if isempty(sym),
         sym = symbol(arg);
      elseif ~isempty(symbol(arg)) & ~strcmp(symbol(arg), sym),
         inconsym  = 1;
         symbol(arg,sym);
      end;
      pn{i} = arg;
   end;
   
   % denominators
   eval('arg = pol(pd{i});', 'error(peel(lasterr));');
   if size(arg) ~= [1 1],
      error('Denominator polynomials must be scalars.');
   end;
   if ~(isempty(symbol(arg)) | strcmp(symbol(arg),'s') | ...
         strcmp(symbol(arg),'p')),
      error('Continuous time only.');
   end;
   ds = max(ds, deg(arg));
   if isempty(sym),
      sym = symbol(arg);
   elseif ~isempty(symbol(arg)) & ~strcmp(symbol(arg), sym),
      inconsym = 1;
      symbol(arg,sym);
   end;
   pd{i} = arg;    
end;

if inconsym
   warning('Polynomials have inconsistent symbols.');
   disp(['Symbol is assumed to be [' sym '].']);
end;

if isempty(d),
 d = ds; % controller order = systems order
elseif ~isa(d,'double') | length(d)~=1 | ...
      ~isreal(d) | d <= 0,
 error('Invalid controller order; must be positive scalar.');
end;

n = 2*d+1; % size of Hermite-Fujiwara matrix
m = ds+d; % degree of characteristic polynomials

if verbose,
 disp(['LMISIMSTAB: Simultaneous stabilization of ' int2str(N) ' systems.']);
 disp(['LMISIMSTAB: Systems order =  ' int2str(ds) '.']);
 disp(['LMISIMSTAB: Controller order = ' int2str(d) '.']);
end;

% Particular cases

% Check if all denominators are stable

stop = 1;
i = 1;
while (i <= N) & stop,
  stop = isstable(pd{i});
  i = i + 1;
end;
if stop, % all denominators are stable: return zero controller
 if verbose,
  disp('LMISIMSTAB: All the systems are already open-loop stable.');
 end;
 nc = pol(0);
 dc = pol(0);
 return;
end;

% Otherwise, build LMI

if verbose,
 disp('LMISIMSTAB: Build polynomials.');
end;

p = {};
for i = 1:N,
 for j = 0:d,
  p{i,j+1} = shift(pn{i}, j);
  p{i,j+d+2} = shift(pd{i}, j);
 end;
end;

if verbose,
 disp('LMISIMSTAB: Build Hermite-Fujiwara matrices.');
end;

H = {};
for i = 1:N,
 for j = 0:n,
  for k = 0:n,
   H{i,j+1,k+1} = hermfuji(p{i,j+1}, p{i,k+1}, m);
  end;
 end;
end;

if verbose,
 disp('LMISIMSTAB: Build LMI.');
end;

setlmis([]);

% vector of components of X
% x = [X11 X21 X22 X31 X32 X33 ... ]
nvar = (n+1)*(n+2)/2;
varx = zeros(nvar, 1);
structX = zeros(n+1);
k = 0;
for i = 1:n+1,
 for j = i:n+1,
  k = k + 1; varx(k) = lmivar(1,[1 1]);
  structX(i,j) = varx(k);
  structX(j,i) = varx(k);
 end;
end;

% stability constraints for each system
for i = 1:N,
 l = 0;
 LMI = newlmi;
 for j = 1:n+1,
  for k = j:n+1,
    l = l + 1;
    if j == k,
     lmiterm([-LMI 1 1 varx(l)], H{i,j,j}, 1);
    else
     lmiterm([-LMI 1 1 varx(l)], H{i,j,k}+H{i,k,j}, 1);
    end; 
  end;
 end;
 lmiterm([LMI 1 1 0], 1e4*tol); % strictly positive definite
end;

% X = symmetric matrix of size n+1, built from components of vector x
varX = lmivar(3, structX);

LMISYS = getlmis;

if verbose,
 disp('LMISIMSTAB: Call LMIRANK to find a rank-one solution to the LMI.');
end;

[X, x] = lmirank(LMISYS, varX, 1, tol, 'random');

if isempty(X),

 if verbose,
  disp('LMISIMSTAB: No rank-one solution was found.');
 end;
 nc = []; dc = [];

else,

 if verbose,
  disp('LMISIMSTAB: A rank-one solution was found.');
 end;

 nc = pol(x(1:d+1)', d);
 dc = pol(x(d+2:2*d+2)', d);

 % Check LMIs and closed-loop stability

 if verbose,

  disp('LMISIMSTAB: Check LMIs with matrix X.');

  for i = 1:N,
   l = 0;
   LMI = zeros(m);
   for j = 1:n+1,
    for k = 1:n+1,
     LMI = LMI + H{i,j,k}*X(j,k);
    end;
   end;
   vpmin = min(eig(LMI));
   if vpmin < 0,
    warning(['LMI #' int2str(i) ' is not satisfied with matrix X !']);
  end;
  end;

  disp('LMISIMSTAB: Check LMIs with vector x.');

  for i = 1:N,
   l = 0;
   LMI = zeros(m);
   for j = 1:n+1,
    for k = 1:n+1,
     LMI = LMI + H{i,j,k}*x(j)*x(k);
    end;
   end;
   vpmin = min(eig(LMI));
   if vpmin < 0,
    warning(['LMI #' int2str(i) ' is not satisfied with vector x !']);
   end;
  end;

  disp('LMISIMSTAB: Check stability of polynomials.');

  for i = 1:N,
   delta = nc*pn{i}+dc*pd{i};
   tH = hermfuji(delta);
   vpmin = min(eig(tH));
   if vpmin < 0,
    warning(['Characteristic polynomial #' int2str(i) ' is unstable !']);
   end;
  end;

 end;
end;

%end .. lmisimstab
