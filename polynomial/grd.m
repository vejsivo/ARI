function [G,U] = grd(varargin)
%GRD  Greatest right divisor
%
% The commmand
%    G = GRD(N1, N2, .., Nk) 
% computes a greatest right common divisor G of several polynomial 
% matrices N1, N2, .., Nk (k > 0) that all have the same number
% of columns.
%
% The rows of G form a polynomial basis for the module spanned by
% the rows of N0 = [N1' N2' .. Nk']'. Note that this basis is not
% necessarily of minimum degree. The number of rows in G is equal to
% the rank of N0. If N0 has full column rank then G is square.
%
% The command
%    GRD(N1, .. , 'gau') 
% computes the divisor through a modified version of Gaussian elimination. 
% This method is preferable esthetically and generally results in a divisor
% of low degree. However, it may not be numerically reliable. This is the 
% default method.
%
% The command
%    GRD(N1, .. , 'syl')
% computes the divisor through stable reductions of Sylvester matrices. This 
% method is preferable numerically but may result in a divisor of high degree.
%
% Matrices Mi such that Ni = Mi*G may be recovered with Mi = XAB(G, Ni).
%
% The commmand
%    [G, U] = GRD(N1, .., Nk) 
% additionally also returns a unimodular matrix U such that 
%    U * N0 = [ G' 0 .. 0 ]'
% For two input matrices N1 and N2, U may be split into two pairs of left 
% coprime polynomial matrices (P,Q) and (R,S) such that
%    U = [ P R |      P*N1 + Q*N2 = G
%        | Q S ]      R*N1 + S*N2 = 0
%
% A tolerance TOL may be specified as an additional input argument.
% To be distinguishable from other arguments, it must have a form of
% a character string, e.e. '1e-8' . Its default value is the global
% zeroing tolerance.
%
% See also: GLD.

%     Author: D. Henrion, March 10, 1999.
%     Updated to 3.0 by D. Henrion, August 30, 2000.
%     Copyright 1998-2000 by Polyx, Ltd.
%     Modified by J.Jezek, 08-Aug-2001, arg checking, empty case
%                 J.Jezek, 19-Oct-2002, nargin
%                 M.Hromcik,9-Aug-2005, modif for M7-R14-SP2

global PGLOBAL;
eval('PGLOBAL.VERBOSE;', 'painit;');

ni = nargin;
if ni==0,
   error('Not enough input arguments.');
end;
n = 0; % number of input matrices
m = 0; % column dimension

% verbose level
verbose = strcmp(PGLOBAL.VERBOSE, 'yes');

% Handle the case when GRD is called through GLD with
% transposed arguments.

if isa(varargin{1}, 'cell'),
 varargin = varargin{1};
end;
ni = length(varargin);

% Parse input arguments.
% Arrange input matrices into a cell array N = {N1, N2, ..}.

N = {}; rN = [];
tol = [];
method = 'gau';
n = 0; degree = -Inf;

for i = 1:ni,
   arg = varargin{i};
   if isa(arg,'char'),
      numarg = str2num(arg);
      if isempty(numarg),
         method = arg;
      else
         tol = numarg;
      end;
   else
      if isa(arg,'double') & ndims(arg)==2,
         arg = pol(arg);
      end;
      if ~isa(arg,'pol'),
         error(['Invalid ',nth(i),' argument.']);
      end;
      if any(any(isnan(arg))) | any(any(isinf(arg))),
         error('Polynomial is not finite.');
      end;
      if n == 0,
         arg1 = arg;
         [p,m] = size(arg); 
         Pv = arg.var; Ph = arg.h;
         w = logical(0); wz = logical(0); wzi = logical(0);
         wh = logical(0);
      else
         [p,q] = size(arg);
         if q ~= m, 
            error('Matrices of inconsistent dimensions.');
         end;
         v = arg.var; h = arg.h;
         if ~isempty(v),
            if isempty(Pv),
               Pv = v;
               if isempty(Ph),
                  Ph = h;
               end;
            elseif ~strcmp(v,Pv),
               if strcmp(v,'z'),
                  if strcmp(Pv,'z^-1') | wzi,
                     error('Inconsistent variables.');
                  else wz = logical(1);
                  end;
               elseif strcmp(v,'z^-1'),
                  if strcmp(Pv,'z') | wz,
                     error('Inconsistent variables.');
                  else wzi = logical(1);
                  end;
               else
                  if strcmp(Pv,'z'), wz = logical(1);
                  elseif strcmp(Pv,'z^-1'), wzi = logical(1);
                  else w = logical(1);
                  end;
               end;
            elseif ~isempty(Ph) & isfinite(Ph) & ...
                  ~isempty(h) & isfinite(h) & h~=Ph,
               wh = logical(1);
            end;
         end;
      end;
      n = n+1;
      rN(n) = p; N{n} = arg;
      degarg = deg(arg);
      if ~isempty(degarg),
         degree = max(degree, degarg);
      end;
   end;
end;

if w | wz | wzi
   warning('Inconsistent variables.');
elseif wh
   warning('Inconsistent sampling periods.');
end;

if isempty(tol),
   tol = PGLOBAL.ZEROING;
elseif length(tol)~=1 | ~isreal(tol) | tol<0 | tol>1,
   error('Invalid tolerance.');
end;

% compute greatest divisor

if strcmp(method, 'syl'),

 % ****************
 % SYLVESTER MATRIX
 % ****************

 % Extract minimum degree polynomial basis via reduction to row staircase
 % form of the compound matrix N0 = [N1; N2; ..]

 if verbose,
  disp('GRD: Divisor extraction by Sylvester matrix approach.');
  disp('GRD: Reduction to staircase form.');
 end;

 srN = cumsum(rN);
 N0 = pol(zeros(srN(n), m));
 for i = 1:n,
  N0(srN(i)-rN(i)+1:srN(i), :) = N{i};
 end;

 [G, U, colind] = tri(N0, 'row', 'syl', tol);
 rk = sum(colind & colind); G = G(1:rk, :);

elseif strcmp(method, 'gau'),

 % *********************
 % GAUSSIAN ELIMINATION
 % *********************

 if verbose,
  disp('GRD: Divisor extraction by Gaussian elimination.');
 end;
 
 p = sum(rN);
 
 if degree < 0
  G = pol(zeros(0,m));
  if nargout > 1,
    U = pol(eye(p));
  end;
  return;
 end;

 % Build Barnett's resultant S.

 q = 0;
 rS = p; cS = (degree+1)*m; S = zeros(rS, cS);

 for i = 1:n,
  rrN = size(N{i}, 1); degN = N{i}.degree;
  if degN >= 0,
   Coef = N{i}.coef;
   S(1+q:rrN+q, 1+(degree-degN)*m:(degree+1)*m) = reshape(Coef(:,:,degN+1:-1:1), size(S(1+q:rrN+q, 1+(degree-degN)*m:(degree+1)*m)));
  end; 
  q = q + rrN;
 end;

 % ------------------------------------------------------------
 % first step : reduction to reduced row echelon form
 % row permutations are allowed
 % ------------------------------------------------------------

 if verbose,
  disp('GRD: Initial reduction to row echelon form.');
 end;

 % tolerance for Gaussian elimination
 me = max(rS, cS) * norm(S, 'inf');
 if isempty(tol),
  tolgauss = me * PGLOBAL.ZEROING;
 else
  tolgauss = me * tol;
 end;

 row = 1; col = 1; colind = [];
 while (row <= rS) & (col <= cS),
  [val,i] = max(abs(S(row:rS,col))); i = i+row-1;
  if (val <= tolgauss), % negligible column
   S(row:rS,col) = zeros(rS-row+1,1);
  else
   colind = [colind col];
   S([row i],col:cS) = S([i row],col:cS); % row permutation
   S(row,col:cS) = S(row,col:cS)/S(row,col); % pivot
   for i = [1:row-1 row+1:rS],
    S(i,col:cS) = S(i,col:cS) - S(i,col)*S(row,col:cS);
   end;
   row = row + 1;
  end;
  col = col + 1;
 end;

 % ------------------------------------------------------------
 % second step : successive reductions without row permutations
 % ------------------------------------------------------------

 E = S; rk = rank(S, tolgauss); oldrk = 0;

 % higher tolerance for Gaussian elimination
 me = max(rS, cS) * norm(S, 'inf') * 1e2;
 if isempty(tol),
  tolgauss = me * PGLOBAL.ZEROING;
 else
  tolgauss = me * tol;
 end;

 iteration = 0;
 first = 0; % enters the loop for the first time
 while ~first | (rk-oldrk > m),

  iteration = iteration + 1;
  if verbose,
   fprintf('GRD: Iteration %3d ', iteration);
  end;

  first = 1;
  oldrk = rk;
  S = [S zeros(rS, m); zeros(p, m) E]; % extended resultant
  rS = rS + p; cS = cS + m; degree= degree + 1;

  % higher tolerance for Gaussian elimination
  me = max(rS, cS) * norm(S, 'inf') * 1e-4;
  if isempty(tol),
   tolgauss = me * PGLOBAL.ZEROING;
  else
   tolgauss = me * tol;
  end;

  % "shifted" row echelon form by Gaussian elimination
  % row permutations are not allowed

  zp = zeros(rS, 1);
  col = 1;
  while (col <= cS), 
   row = 1;
   while (row <= rS),
    if (abs(S(row,col)) > tolgauss) & ~zp(row), % non-zero possible pivot
     zp(row) = 1; % next pivot shouldn't belong to this row
     S(row,:) = S(row,:)/S(row,col); % normalization
     for k = 1:rS, % cancel other elements in the column
      if (k ~= row) & (abs(S(k,col)) > tolgauss),
       S(k,:) = S(k,:) - S(row,:)*S(k,col);
      end;
     end;
     row = rS;
    end; % if
    row = row + 1;
   end; % while row
   col = col + 1;
  end; % while column

  rk = rank(S, tolgauss);
  E = S(rS-p+1:rS,:); % lower part of S

 if verbose,
  fprintf('   Number of missing rows = %3d.\n', rk-oldrk-m);
 end;

 end; % while rk

 % --------------
 % grd extraction
 % --------------

 ind = any(abs(E')>tolgauss);
 if sum(ind),
  E = E(ind, :); % non zero rows
  for i=0:degree,
   G(:, 1+i*m:(i+1)*m) = E(:, 1+cS-m*(i+1):cS-m*i);
  end;  
 else
  error('Extraction failed. Try to modify tolerance.');
 end;

 % zero leading coefficients will be removed
 G = pzer(pol(G, degree), tolgauss);
 G.v = Pv; G.h = Ph;

 % check whether G is square
 rowG = size(G, 1); rankG = rank(G, tolgauss);
 if (rowG ~= m) | (rankG < m),

  % extraction of the linearly independent
  % rankG rows of G(s) with minimal row degrees

  [rowdeg, ind] = sort(deg(G, 'row'));
  newG = pol(zeros(0,size(G,2)));
  oldrank = 0; j = 0;
  for i = 1:rankG,
   newrank = oldrank;
   while (newrank == oldrank) & (j < rowG),
    j = j + 1;
    testG = [newG; G(ind(j), :)];
    newrank = rank(testG, tolgauss);
   end;
   if newrank > oldrank,   
    newG = testG;
    oldrank = newrank;
   end;
  end;
  G = newG;
 end;

 % ------------------------
 % Recover reduction matrix
 % ------------------------

 if nargout > 1,

  if verbose,
   disp('GRD: Recover reduction matrix.');
  end;

  % First rows of reduction matrix are recovered by solving
  % a system of equations: U1*[N1;N2..] = G

  srN = cumsum(rN);
  N0 = pol(zeros(srN(n), m));
  for i = 1:n,
     arg = N{i};
     arg.v = v; arg.h = h;
     N0(srN(i)-rN(i)+1:srN(i), :) = arg;
  end;

  U1 = xab(N0, G, tol);
  
  if any(any(isnan(U1))),
   disp('GRD: Failed to recover reduction matrix.');
   error('Try to modify tolerance.');
  end;

  % Remaining rows are recovered by extracting the left
  % null-space of matrix [N1;N2..], i.e. U2*[N1;N2..] = 0

  U2 = null(N0.', tol).';

  if any(any(isnan(U2))),
   disp('GRD: Failed to recover reduction matrix.');
   error('Try to modify tolerance.');
  end;

  U = vertcat(U1,U2);
   
 end;

else

 error('Invalid command option.');

end;

if verbose,
 disp(['GRD: A divisor of degree ' int2str(deg(G)) ' was found.']);
end;

%end .. grd

