function [a,b,c,d] = rmf2ss(N,D,tol)
%RMF2SS  Convert a right matrix fraction to a controller-form realization
%
% Given two polynomial matrices N and D such that D is square nonsingular
% column reduced the command
%    [a,b,c,d] = RMF2SS(N, D [,TOL])
% returns the (generalized) controller-form realization (a,b,c,d) of the
% system with transfer matrix H, that is,
%     H(x) = N(x) D^-1(x) = c (xI-a)^-1 b +d(x),
% where x is 's' (continuous-time) or 'z' (discrete-time).
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.
%
% See also: LMF2SS.

%    Author: D. Henrion, M. Sebek, R.C.W. Strijbos, December 11, 1998.
%    $ Revision 3.0 $   $ Date 08-Jun-2000  J.Jezek  $%
%    Last modified by D. Henrion, August 29, 2000.
%                  by J. Jezek, Jul 22, 2001, Var, sampl per check $
%    Copyright 1998-2000 by Polyx, Ltd.
    
global PGLOBAL;
eval('PGLOBAL.ZEROING;', 'painit;');

% verbose level
verbose = strcmp(PGLOBAL.VERBOSE, 'yes');

if nargin < 2
   error('Not enough input arguments.');
end
eval('N = pol(N); D = pol(D);', ...
   'error(peel(lasterr));');

eval('[N,D] = testdnd(N,D,''r'');', ...
   'error(peel(lasterr));');

if nargin > 2 & ~isempty(tol),
   if ~isa(tol, 'double') | (length(tol) > 1) | ...
         ~isreal(tol) | tol<0 | tol>1,
      error('Invalid tolerance.');
   end;
else
   tol = PGLOBAL.ZEROING; % default tolerance
end;

[rN, cN] = size(N);
[rD, cD] = size(D);
% Empty matrices: section added by D. Henrion, August 1, 2000
if (rN == 0) | (cN == 0), % empty matrices
 a = zeros(0); b = zeros(0,cN); c = zeros(rN,0); d = zeros(rN,cN);
 return;
end;

Var = ''; h = 0;
eval('[Var,h,N,D] = testvhnd(N,D);', ...
   'error(peel(lasterr));');

if strcmp(Var, 'z^-1') | strcmp(Var, 'd'),
   [N,D] = reverse(N,D,'r',tol);
   N.v = 'z'; D.v = 'z';
end;

dcolN = deg(N, 'col'); dcolN(dcolN < 0) = 0;
dcolD = deg(D, 'col'); dcolD(dcolD < 0) = 0; % controllability indices
colD = lcoef(D, 'col');

d = pol(zeros(rN,rD));
S = svd(colD);
if min(S) < rD * max(S) * tol,
   error('Denominator matrix is singular or not column reduced.');
elseif max(dcolN - dcolD) >=0 %& dcolN ~= 0,
   [d,N]=rdiv(N,D);         % Extraction of the proper part 
end;

n = sum(dcolD); p = rD; m = rN;
a = zeros(n); b = zeros(n,p); c = zeros(m,n);
% computation of block indices in A
ind = [1 dcolD]; ind = cumsum(ind);
invDhc = inv(colD);

% lower degree coefficient matrices Dlc and C = Nlc
Dlc = zeros(p, n);
for j = 1:p,
   for i = 1:p,
      entry = D(i,j); dge = deg(entry);
      if dge < 0,
	 poly = 0;
      else
	 poly = zeros(1, dge+1); poly(:) = entry.coef;
      end;
      poly = [poly zeros(1,dcolD(j)+1-length(poly))];
      Dlc(i,ind(j):ind(j+1)-1) = fliplr(poly(1:length(poly)-1));
   end;
   for i = 1:m,
      entry = N(i,j); dge = deg(entry);
      if dge < 0,
	 poly = 0;
      else
	 poly = zeros(1, dge+1); poly(:) = entry.coef;
      end;
      poly = [poly zeros(1,dcolD(j)+1-length(poly))];
      c(i,ind(j):ind(j+1)-1) = fliplr(poly(1:length(poly)-1));
   end;
end;

% a and b
K = invDhc*Dlc;
for j = 1:p,
   if dcolD(j) > 1,
      a(ind(j):ind(j+1)-1,ind(j):ind(j+1)-1) = diag(ones(1,dcolD(j)-1),-1);
   end;
   if dcolD(j) > 0,
      a(ind(j),:) = - K(j,:);
      b(ind(j),:) = invDhc(j,:);
   end;
end;

%d
if ~isnumeric(d)
   d = pzer(d);
   if d.deg == 0 | isinf(d.deg),
      d=d{:};
   elseif strcmp(Var,'z^-1')
      d.var = 'z';
   elseif strcmp(Var,'d')
      d.var = 'q';
   end
end

%end .. rmf2ss


