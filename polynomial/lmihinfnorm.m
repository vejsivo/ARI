function [gammainf,gammasup] = lmihinfnorm(N,D,opt1,opt2,tol)
%LMIHINFNORM  H-infinity norm of a polynomial matrix fraction via LMIs
%
% The commands
%    
%    GAMMAUP = LMIHINFNORM(N,D)
%    GAMMAUP = LMIHINFNORM(N,D,'l')
%
% compute a tight upper bound on the H-infinity norm of the stable
% rational matrix G = D \ N. Similarly, the command
%
%    GAMMAUP = LMIHINFNORM(N,D,'r')
%
% computes the same thing for the stable rational matrix G = N / D.
% The upper bound is computed via a dual LMI optimization problem.
%
% With two output arguments, the commands
%
%    [GAMMALOW,GAMMAUP] = LMIHINFNORM(N,D)
%    [GAMMALOW,GAMMAUP] = LMIHINFNORM(N,D,'bin')
%
% computes also computes a tight lower bound. The lower bound is computed
% via binary search over successive LMI feasibility problems. With the syntax
%
%    [GAMMALOW,GAMMAUP] = LMIHINFNORM(N,D,'gevp')
%
% the lower bound is computed via a generalized eigenvalue LMI problem.
% The latter option is generally faster than the former, but the generalized
% eigenvalue solver sometimes features slow convergence and poor accuracy.
%
% The relative accuracy for solving LMIs can be specified as an optional
% input argument. Its default value is the global zeroing tolerance.

% Author: Didier Henrion, February 23, 2000.
% Updated to 3.0 by Didier Henrion, September 1, 2000.
%                by Jan Jezek, Aug 2001, arg checking
% Copyright 2000 Polyx, Ltd.

global PGLOBAL;
eval('PGLOBAL.ZEROING;', 'painit;');

if nargin < 2,
   error('Not enough input arguments.');
elseif nargin < 3,
  opt1 = []; opt2 = []; tol = [];
elseif nargin < 4,
  opt2 = []; tol = [];
  if isa(opt1, 'double'), tol = opt1; opt1 = []; end;
elseif nargin < 5,
  tol = [];
  if isa(opt1, 'double'), tol = opt1; opt1 = []; end;
  if isa(opt2, 'double'), tol = opt2; opt2 = []; end;
else
  if isa(opt1, 'double'), swap = tol; tol = opt1; opt1 = swap; end;
  if isa(opt2, 'double'), swap = tol; tol = opt2; opt2 = swap; end;
end;

eval('N = pol(N); D = pol(D);', 'error(peel(lasterr));');

if ~strcmp(opt1, 'r') & ~strcmp(opt2, 'r'),
  N = N.'; D = D.'; % left MFD is transposed
end;

if strcmp(opt1, 'gevp') | strcmp(opt2, 'gevp'),
 binary = 0; % GEVP for lower bound
else
 binary = 1; % binary search for lower bound
end;

if isempty(tol),
   tol = PGLOBAL.ZEROING;
elseif length(tol)~=1 | ~isreal(tol) | tol<0 | tol>1
   error('Invalid tolerance.');
end;

verbose = strcmp(PGLOBAL.VERBOSE, 'yes');

[rN,cN] = size(N);
[rD,cD] = size(D);

if any(any(isnan(N))) | any(any(isnan(D))) | ...
   any(any(isinf(N))) | any(any(isinf(D))),
      error('Polynomial is not finite.');
end;

if rD ~= cD,
 error('Denominator matrix is not square.');
elseif cN ~= rD,
   if (rD == 1),
      D = diag(repmat(D,cN,1));
      rD = cN; cD = cN;
   else
      error('Matrices of inconsistent dimensions.');
   end;
end;

if rN == 0 | cN == 0,   % empty matrices
   gammainf = 0; gammasup = 0;
   return;
end;

[tv,sym,N,D] = testvp(N,D);
if tv==2,
   [th,h,N,D] = testhp(N,D,sym);
   if ~th,
      warning('Inconsistent sampling periods.');
   end;
   if strcmp(D.v,'z^-1'),
      dg = deg(D);
      N = shift(N,dg); D = rev(D,dg);
      D.v = 'z';
   else
      dg = deg(N);
      N = rev(N,dg); D = shift(D,dg);
      N.v = 'z';
   end;
   sym = 'z';
elseif ~tv,
   warning('Inconsistent variables.');
end;

[th,h,N,D] = testhp(N,D,sym);
if ~th,
   warning('Inconsistent sampling periods.');
end;

if issingular(D,tol),
   warning('Denominator matrix is singular to working precision.');
   gammainf = Inf;
   gammasup = Inf;
   return;
end;

dN = deg(N);
dD = deg(D);

% number of digits to be displayed

ndigits = max(ceil(-log10(tol))+1,3);

% stability region

switch sym,
case {'s','p',''},
   H = [0 1;1 0];
   if dN > dD,
      warning('Fraction is not proper.');
      gammainf = Inf;
      gammasup = Inf;
      return;
   end;
case {'z^-1','d'},
   H = [-1 0;0 1];
case {'z','q'},
   H = [1 0;0 -1];
end;

% projection matrices

PP = [eye(rD*dD) zeros(rD*dD,rD)];
QQ = [zeros(rD*dD,rD) eye(rD*dD)];

% matrix polynomial coefficients

DD = D{:};
NN = [N{:} zeros(rN, rD*(dD-dN))];

% Schur decomposition for trace function

[UD, TD] = schur(DD'*DD);
[UN, TN] = schur(NN'*NN);

% ------------------------------------------
% Upper bound: dual LMI optimization problem
% ------------------------------------------

if verbose,
 disp(['LMIHINFNORM: Relative accuracy = ' num2str(tol)]);
 disp('LMIHINFNORM: Compute upper bound.');
end;

setlmis([]);

vargamma = lmivar(1, [1 1]);
varP = lmivar(1, [rD*dD 1]);

LMI = newlmi;
lmiterm([-LMI 1 1 varP], -H(1,1)*QQ', QQ);
lmiterm([-LMI 1 1 varP], -H(2,1)*PP', QQ, 's');
lmiterm([-LMI 1 1 varP], -H(2,2)*PP', PP);
lmiterm([-LMI 1 1 0], -NN'*NN);
lmiterm([-LMI 1 1 vargamma], DD', DD);

lmiterm([-LMI 1 1 0], eps);
% this small perturbation is really necessary
% otherwise the LMI is sometimes found infeasible !

LMISYS = getlmis;

nd = decnbr(LMISYS);
c = zeros(1,nd); c(1) = 1;

if verbose,
 disp('LMIHINFNORM: Solve dual LMI optimization problem...');
end;

options(1) = tol; % relative accuracy
options(2) = 100; % max number of iterations
options(4) = 20; % stopping criterion
options(5) = 1; % trace off

[copt,xopt] = mincx(LMISYS,c,options);

if isempty(xopt) | (copt < 0),
 disp('LMIHINFNORM: Infeasible dual LMI !');
 disp('LMIHINFNORM: Fraction may have zeros on the stability boundary.');
 gammainf = Inf;
 gammasup = Inf;
 return;
else
 gammasup = sqrt(copt);
end;

if verbose,
 disp(['LMIHINFNORM: Computed upper bound = ' ...
       num2str(gammasup,['%.' int2str(ndigits) 'f']) '.']);
end;

if nargout < 2,

 gammainf = gammasup;

else,

 % --------------------------------
 % Lower bound: primal GEVP problem
 % --------------------------------

 if verbose,
  disp('LMIHINFNORM: Compute lower bound.');
 end;

 setlmis([]);

 varX = lmivar(1, [rD*(dD+1) 1]);

 % stability constraint

 LMI = newlmi;
 lmiterm([-LMI 1 1 varX], H(1,1)*QQ, QQ');
 lmiterm([-LMI 1 1 varX], H(2,1)*PP, QQ', 's');
 lmiterm([-LMI 1 1 varX], H(2,2)*PP, PP');

 % X >= 0
 LMI = newlmi;
 lmiterm([-LMI 1 1 varX], 1, 1);

 % found feasible point to initialize GEVP problem

 LMISYS = getlmis;

 setlmis(LMISYS);

 % compute lower found gammainf via SVD of the TF matrix performed
 % at complex points on the stability boundary that are functions of
 % the zeros of the denominator of the MFD
 gammainf = 0;

 if verbose,
  disp('LMIHINFNORM: Compute initial lower bound using SVD.');
 end;
  
 rootD = roots(D, 'eig');

 if H(1,1) == 0, % ct systems
  rootD = [0; sqrt(-1)*abs(rootD)]; % static + SVD at j*abs(zero)
 else % dt systems
  rootD = [1; exp(sqrt(-1)*angle(rootD))]; % static + SVD at exp(j*angle(zero))
 end;

 for i = 1:length(rootD),  
 
  Dv = polyval(D, rootD(i));
  Nv = polyval(N, rootD(i));

  if rank(Dv) < rD,
   disp('LMIHINFNORM: Fraction has zeros on the stability boundary.');
   disp('LMIHINFNORM: H-inf norm is infinite.');
   gammainf = Inf;
   gammasup = Inf;
  else
   gammainf = max(gammainf, max(svd(Nv/Dv)));
  end;

 end;

 if verbose,
  disp(['LMIHINFNORM: Computed initial lower bound = ' ...
       num2str(gammainf,['%.' int2str(ndigits) 'f']) '.']);
 end;

 if gammainf > gammasup,
  error('Lower bound is greater than upper bound !');
 end;

 LMISYS = getlmis;

 if binary,

  % ------------------------------------------------------------
  % Solve GEVP via a binary search between gammainf and gammasup

  if verbose,
   disp('LMIHINFNORM: Perform binary search..');
  end;

  options(1) = tol; % relative accuracy
  options(2) = 100; % max number of iterations
  options(4) = 10; % stopping criterion
  options(5) = 1; % trace off

  while (gammasup-gammainf)/gammasup > tol,

   if verbose,
    fprintf(['LMIHINFNORM: Bounds = [%.' int2str(ndigits) ...
             'f,%.' int2str(ndigits) 'f].\n'], gammainf, gammasup);
   end;

   gammanew = (gammainf+gammasup)/2;

   setlmis(LMISYS);

   % new constraint Trace N'*N*X >= gammainf^2 * Trace D'*D*X

   LMI = newlmi;

   for i = 1:(dD+1)*rD,
    lmiterm([-LMI 1 1 varX], TN(i,i)*UN(:,i)', UN(:,i));
    lmiterm([-LMI 1 1 varX], -gammanew^2*TD(i,i)*UD(:,i)', UD(:,i)); 
   end;
 
   LMISYS_FEAS = getlmis;

   [tmin,x0] = feasp(LMISYS_FEAS,options);

   if isempty(x0),
    gammasup = gammanew;
   else
    gammainf = gammanew;
   end;

  end; % while

 else,

   % -------------------------------------
   % Solve GEVP via LMI toolbox macro GEVP

   if verbose,
    disp('LMIHINFNORM: Compute initial feasible point for GEVP...');
   end;

   setlmis(LMISYS);
   
   % new constraint Trace N'*N*X >= gammainf^2 * Trace D'*D*X

   % slightly less than gammainf to allow warm start in GEVP
   gammagevp = 0.99*gammainf;

   LMI = newlmi;

   for i = 1:(dD+1)*rD,
    lmiterm([-LMI 1 1 varX], TN(i,i)*UN(:,i)', UN(:,i));
    lmiterm([-LMI 1 1 varX], -gammagevp^2*TD(i,i)*UD(:,i)', UD(:,i)); 
   end;
 
   LMISYS_FEAS = getlmis;

   % compute initial feasible point for gammainf

   options(1) = tol; % relative accuracy
   options(2) = 100; % max number of iterations
   options(4) = 10; % stopping criterion
   options(5) = 1; % trace off

   [tmin,x0] = feasp(LMISYS_FEAS,options);

   if isempty(x0),
    disp('LMIHINFNORM: Unfeasible primal LMI problem !');
    disp('LMIHINFNORM: Fraction may have zeros on the stability boundary.');
    gammainf = Inf;
    gammasup = Inf;
    return;
   end;

   if verbose,
    disp('LMIHINFNORM: Solve primal GEVP problem...');
   end;

   setlmis(LMISYS);

   % Trace N'*N*X >= gamma^2 * trace D'*D*X
   LMI = newlmi;

   for i = 1:(dD+1)*rD,
    lmiterm([LMI 1 1 varX], -TN(i,i)*UN(:,i)', UN(:,i));
    lmiterm([-LMI 1 1 varX], TD(i,i)*UD(:,i)', UD(:,i)); 
   end;

   LMISYS = getlmis;

   options(1) = tol; % relative accuracy
   options(2) = 500; % max number of iterations
   options(4) = 10; % stopping criterion
   options(5) = 1; % trace off
   target = -1e20;

   [tmin,xopt] = gevp(LMISYS,1,options,-gammagevp^2,x0,target);

   if tmin > 0,
    disp('LMIHINFNORM: Negative optimal value of primal LMI problem !');
    disp('LMIHINFNORM: Fraction may have zeros on the stability boundary.');
    gammainf = Inf;
    gammasup = Inf;
    return;
   end;

   gammainf = sqrt(-tmin);

 end; % binary search or GEVP problem

 if verbose,
  disp(['LMIHINFNORM: Computed lower bound = ' ...
       num2str(gammainf,['%.' int2str(ndigits) 'f']) '.']);
 end;

end; % if nargout > 1
 
%end .. lmihinfnorm
 

