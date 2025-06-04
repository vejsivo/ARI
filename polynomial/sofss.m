function inter = sofss(q,p,S)
%SOFSS Static output feedback simultaneous stabilization of scalar plants
%
% Given 2 cell arrays of polynomials A = {A1...AN} and B = {B1...BN}, the
% instruction
%  
%   INTER = SOFSS(A,B)
%
% computes the open intervals INTER = [[Kinf1 Ksup1];..;[KinfM KsupM]] such
% that any scalar K chosen within the intervals is a stabilizing static
% output feedback gain, i.e. such that every characteristic polynomial
% Ai + K*Bi is stable
%
% With the optional input argument
%
%   INTER = SOFSS(A,B,S)
%
% one can specify an Hermitian 2x2 stability region matrix S such that
% each closed-loop pole Z must satisfy the quadratic inequality
% S(1,1)+S(1,2)*Z+S(2,1)*Z'+S(2,2)*Z*Z' < 0. Standard choices are
% S=[0 1;1 0] (left half-plane) and S=[-1 0;0 1] (unit disk).
% The default value of S is determined by the variable symbol of input
% polynomials in cell arrays A and B.

% Author: D. Henrion, January 12, 2001
% Last modified by D. Henrion, October 29, 2002
% Copyright 2001-02 by PolyX, Ltd.
  
% Parse input arguments
  
global PGLOBAL;

eval('PGLOBAL.FORMAT;',...
     'error(''Use PINIT to initialize the Polynomial Toolbox.'');');

if nargin < 2,
 error('Invalid number of input arguments');
end;

if ~isa(p, 'cell'), p = {p}; end;
if ~isa(q, 'cell'), q = {q}; end;

tol = PGLOBAL.ZEROING;

verbose = strcmp(PGLOBAL.VERBOSE, 'yes');

N = length(p);
if (nargin > 1),
 if length(q) ~= N,
  error('Invalid input arguments');
 end;
end;

% degree
n = -Inf;
for i = 1:N,
 p{i} = pol(p{i});
 q{i} = pol(q{i});
 n = max([n deg(p{i}) deg(q{i})]);
end;

% variable symbol
symb = [];
for i = 1:N,
 sp = symbol(p{i}); sq = symbol(q{i});
 if isempty(symb), symb = sp; end;
 if isempty(symb), symb = sq; end;
 if ~isempty(symb),
  if (~isempty(sp) & ~strcmp(sp, symb)) | ...
     (~isempty(sq) & ~strcmp(sq, symb)),
   error('Inconsistent variable symbolic in input polynomials');
  end;
 end;
end;

% stability region
if nargin < 3,
 S = [];
end;
if ~isempty(S),
 if ~isa(S,'double'),
   error('Invalid third input argument.');
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
 
% Build quadratic matrix H{i} for each plant

if verbose,
  disp(['SOFSS: ' int2str(N) ' plants of order ' int2str(n) '.']);
  disp('SOFSS: Build quadratic Hermite matrices.');
end;

H = cell(N);
C = hermcoef(S,n);
for i = 1:N,
  H{i} = pol(zeros(n));
  j = 0;
  for i1 = 0:n,
    for i2 = 0:i1,
      j = j+1;
      HC = trimat(C(:,j));
      if norm(HC,1) > 0,
       H{i} = H{i} + (q{i}{i1}+s*p{i}{i1})*(q{i}{i2}+s*p{i}{i2})*HC;
      end;
    end;
  end;
end;

% Extract the real roots of each H{i}

if verbose,
  disp('SOFSS: Extract and sort the real roots of the Hermite matrices.');
end;

r = [];
for i = 1:N
  r = [r; roots(H{i})];
end;
r = real(r(abs(imag(r)) < 1e3*tol)); % purely real zeros only

% Sort the distinct real zeros
r = sort([-Inf; r; +Inf]);
%r = sort([-Inf; r; 0; +Inf]);
i = 2;
while i <= length(r),
 if abs(r(i-1)-r(i)) < 1e3*tol,
   r = r([1:i-1 i+1:length(r)]);
 else
   i = i+1;
 end;
end;

% Evaluate the inertia in each interval

if verbose,
  disp('SOFSS: Evaluate the inertia in each interval.');
end;

npos = zeros(length(r));
for i = 1:length(r)-1,
  if isinf(r(i)), % -Inf
    if isinf(r(i+1)); % +Inf
      x = 0;
    elseif r(i+1) < 0,
      x = 2*r(i+1);
    else
      x = 0;
    end;
  elseif isinf(r(i+1)), % +Inf
    if r(i) < 0,
      x = 0;
    else
      x = 2*r(i);
    end;
  else
    x = (r(i)+r(i+1))/2;
  end;
  rx = [];
  for j = 1:N,
    rx = [rx; real(eig(polyval(H{j}, x)))];
  end;
  npos(i) = sum(rx>tol); % number of positive eigenvalues

  if verbose,
   disp(['In [' num2str(r(i)) ',' num2str(r(i+1)) '] ' num2str(npos(i)) ...
	' positive eigenvalues over ' num2str(n*N) '.']);
  end;

end;

% Extract the intervals where eigenvalues are all positive
inter = []; j = 0;
if npos(1) == n*N,
  j = j+1; inter(j,:) = [r(1) r(2)];
end;
for i = 2:length(r)-1,
  if npos(i) == n*N,
    if npos(i-1) ~= n*N,
     j = j+1; inter(j,:) = [r(i) r(i+1)];
    else
     inter(j,2) = r(i+1);
    end;
  end
end;

% END OF FUNCTION SOFSS.M

function mat = trimat(vec)
m = length(vec);
n = (1+sqrt(1+8*m))/2-1;
mat = zeros(n);
row = 1;
for rowtri = 1:n,
  mat(rowtri, 1:rowtri) = vec(row:row+rowtri-1)';
  row = row+rowtri;
end;
mat = mat+mat'-diag(diag(mat));

