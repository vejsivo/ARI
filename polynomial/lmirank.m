function [X,Y] = lmirank(LMISYS, varX, rankX, tol, opt)
%LMIRANK  Solve an LMI for a matrix of specified rank
%
% The command
%
%   [X,Y] = LMIRANK(LMISYS, VARX, RANKX, TOL)
%
% seeks a symmetric positive semi-definite matrix X of rank RANKX
% solution to a system of LMIs built using the LMI Control Toolbox.
% The system of LMIs is specified through the internal description
% LMISYS returned by GETLMIS. The identifier of matrix X within this
% system is specified through VARX.
%
% The macro also returns a full-rank matrix Y such that X = Y*Y'. 
%
% Solving an LMI with rank constraint is a difficult non-convex problem.
% This macro is just an attempt to solve this problem via a heuristic
% based on convex relaxations, more specifically an implementation of
% a cone-complementarity algorithm. The heuristic is not guaranteed to
% provide a solution, even when one exists. If the macro fails to find a
% solution, then X is empty.
%
% If S denotes the vector of singular values of matrix X (as given by
% the command S = SVD(X)), then the algorithm stops
% - either if S(RANKX) / S(RANKX+1) is less than relative accuracy TOL,
%   in which case a solution X of RANKX is found,
% - or if the number of iterations is greater than 15, in which case no
%   solution X of RANKX is found.
% The default value of TOL is the global zeroing tolerance.

% The cone-complementarity algorithm is described in L. El Ghaoui, F. Oustry,
% M. Ait Rami "A Cone Complementariy Linearization Algorithm for Static
% Output Feedback and Related Problems" IEEE Trans. Autom. Control, Vol. 42,
% No. 8, pp. 1171-1176, 1997. Its application to the rank minimisation
% under LMI constraints is described e.g. in Chapter 1 in L. El Ghaoui,
% S. Niculescu (Editors) "Advances in LMI Methods in Control", SIAM,
% Philadelphia, 1999.

% Author: D. Henrion, January 25, 2000.
% Modified by D. Henrion, February 11, 2000.
% Updated to 3.0 by D. Henrion, September 1, 2000.
% Copyright 2000 by Polyx, Ltd.

global PGLOBAL
eval('PGLOBAL.ZEROING;', 'painit;');

NBMAXITER = 15; % total number of iterations

if nargin < 3,
 error('Not enough input arguments.');
elseif nargin < 4,
 tol = []; opt = [];
elseif nargin < 5,
 opt = [];
end;

if isa(tol, 'char'),
 swap = tol; tol = opt; opt = swap;
end;

if isempty(tol),
 tol = PGLOBAL.ZEROING;
end;

if ~strcmp(opt, 'random'),
 opt = [];
end;

verbose = strcmp(PGLOBAL.VERBOSE, 'yes');

% ---------------------------------
% Add constraints to the LMI system
% ---------------------------------

setlmis(LMISYS);

X = decinfo(LMISYS, varX);
n = size(X, 1);

if rankX > n,
 error('The expected rank exceeds the dimension of the matrix variable.');
elseif rankX < 1,
 error('The expected rank must be a strictly positive integer.');
end;

varvx = lmivar(2, [n rankX]);

varY = lmivar(1, [n 1]);
varvy = lmivar(2, [n rankX]);

varZ = lmivar(1, [rankX 1]);

% constraint [X x;x' 1] >= 0
LMI = newlmi;
lmiterm([-LMI 1 1 varX], 1, 1);
lmiterm([-LMI 1 2 varvx], 1, 1);
lmiterm([-LMI 2 2 0], 1);

% constraint [Y y;y' z] >= 0
LMI = newlmi;
lmiterm([-LMI 1 1 varY], 1, 1);
lmiterm([-LMI 1 2 varvy], 1, 1);
lmiterm([-LMI 2 2 varZ], 1, 1);

% constraint Y >= I
LMI = newlmi;
lmiterm([-LMI 1 1 varY], 1, 1);
lmiterm([LMI 1 1 0], 1);

LMISYS = getlmis;

nvar = decnbr(LMISYS);

options(1) = tol; % relative accuracy on criterion
options(5) = 1; % trace off

% ---------------------------------
% Iterative Frank and Wolfe process
% ---------------------------------

% First step: feasible point X0, Y0

if isempty(opt),

 if verbose,
  disp('LMIRANK: Compute initial point.');
 end;

 % Solve LMI feasibility problem
 [tmin, xfeas] = feasp(LMISYS, options);

 if isempty(tmin) | tmin > 0,
  error('LMI problem is infeasible !');
 end;

else

 % Random start: solve LMI optimization problem with random criterion

 if verbose,
  disp('LMIRANK: Compute a random initial point.');
 end;

 c = randn(decnbr(LMISYS), 1);
 [copt, xfeas] = mincx(LMISYS, c, options);

 if isempty(xfeas),
  error('LMI problem is infeasible !');
 end;

end;

% Retrieve variables
X = dec2mat(LMISYS, xfeas, varX); vx = dec2mat(LMISYS, xfeas, varvx);
Y = dec2mat(LMISYS, xfeas, varY); vy = dec2mat(LMISYS, xfeas, varvy);
Z = dec2mat(LMISYS, xfeas, varZ);
XX = [X vx;vx' eye(rankX)];
YY = [Y vy;vy' Z];

% Next steps: minimization of a linearized criterion

nit = 0;
S = [svd(X); 0]; criterion = S(rankX+1) / S(rankX);

while (criterion > tol) & (nit < NBMAXITER),

 nit = nit + 1;

 if verbose,
  fprintf('LMIRANK: Iteration %2d ..', nit);
 end;

 % Linearized criterion for LMI = trace(XX*YYk+YY*XXk)
 c = zeros(nvar, 1);
 for i = 1:nvar,
  [Xi, vxi, Yi, vyi, Zi] = defcx(LMISYS, i, varX, varvx, varY, varvy, varZ);
  XXi = [Xi vxi; vxi' eye(rankX)];
  YYi = [Yi vyi; vyi' Zi];
  c(i) = trace(XXi*YY) + trace(YYi*XX);
 end;

 % Solve LMI optimization problem

 [copt, xopt] = mincx(LMISYS, c, options);

 if isempty(copt),

  disp('WARNING: LMI problem becomes infeasible !');
  nit = NBMAXITER;

 else

  % Retrieve matrix variables
  X = dec2mat(LMISYS, xopt, varX); vx = dec2mat(LMISYS, xopt, varvx);
  Y = dec2mat(LMISYS, xopt, varY); vy = dec2mat(LMISYS, xopt, varvy);
  Z = dec2mat(LMISYS, xopt, varZ);
  XX = [X vx;vx' eye(rankX)]; YY = [Y vy;vy' Z];

  % Stopping criterion: ratio between singular value #rankX+1 and
  % singular value #rankX must be small enough

  S = [svd(X); 0]; criterion = S(rankX+1) / S(rankX);

  if verbose,
   disp([' Criterion = ' num2str(criterion)]);
  end;

 end;

end;

Y = vx;

if nit == NBMAXITER,
 disp(['LMIRANK: The heuristic failed. However, it does not necessarily']);
 disp(['         mean that there is no solution of rank ' int2str(rankX) '.']);
% X = []; Y = [];
Y=vx;
else
 Y = vx;
end;

%end .. lmirank
