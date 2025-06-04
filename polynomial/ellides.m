function [x,y] = ellides(a0,b0,a1,b1,D,H,pid)
% ELLIDES - Robust stabilization of an ellipsoid of polynomials
%  
% Given polynomials A0, B0 and A1, B1 and C the instruction
%
%  [X,Y] = ELLIDES(A0,B0,A1,B1,C)
%
% attempts to find two polynomials X and Y such that the proper
% controller Y/X robustly stabilizes the plant with transfer function
% (B0+D'*B1)/(A0+D'*A1) where D is a real vector such that NORM(D) <= 1
%
% Roughly speaking, polynomial C corresponds to the central (or nominal)
% closed-loop polynomial, for example obtained by stabilizing nominal
% plant B0/A0 (see macro AXBYC). It fixes the order of the sought
% controller.
%
% With the optional input arguments
%
%  [X,Y] = ELLIDES(A0,B0,A1,B1,C,S,PID)
%
% one can specify an Hermitian 2x2 stability region matrix S such that
% each closed-loop pole Z must satisfy the quadratic inequality
% S(1,1)+S(1,2)*Z+S(2,1)*Z'+S(2,2)*Z*Z' < 0. Standard choices are
% S=[0 1;1 0] (left half-plane) and S=[-1 0;0 1] (unit disk).
% The default value of S is determined by the variable symbol of input
% polynomials A0,B0,A1,B1,C.
%
% When the last input argument PID (could be any declared variable or
% constant) is present, then the controller denominator is set to X = s.
% So if central polynomial C has degree 1, the sought controller is
% proportional-integral Y/X = Kp+Ki/s, and if C has degree 2, the
% controller is a PID Y/X = Kp+Ki/s+Kd*s.
%
% The function is based on the theory described in [M. Kvasnica,
% D. Henrion. Low-order robust controller sythesis for systems affected
% by ellipsoidal uncertainty: a polynomial approach, 2002]
%  
% The LMI optimization problem is solved with SeDuMi and its
% LMI interface. Both packages must be properly installed.
 
% Written by Michal Kvasnica and Didier Henrion, February 11, 2002
% Last modified by Didier Henrion, September 16, 2002
% Copyright 2002 by PolyX, Ltd.  

global PGLOBAL;

eval('PGLOBAL.FORMAT;',...
     'error(''Use PINIT to initialize the Polynomial Toolbox.'');');

if nargin < 5,
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

a0 = pol(a0); b0 = pol(b0);
a1 = pol(a1); b1 = pol(b1);
D = pol(D);

n = length(a1)+1;

% d - degree of the denominator, 
d = deg(D);

% variable symbol
symb = symbol(a0);
M = {b0,a1,b1,D};
for i = 1:4,
 newsymb = symbol(M{i});
 if isempty(symb),
  symb = newsymb;
 elseif ~isempty(newsymb) & ~strcmp(symb, newsymb),
  error('Inconsistent variable symbols in input arguments');
 end;
end;

% Check degrees
if deg(b0) > deg(a0),
 error('System is not proper');
end;
if deg(b1) > deg(b0),
 error('Invalid degree of uncertain numerator polynomial');
end;
if deg(a1) > deg(a0),
 error('Invalid degree of uncertain denominator polynomial');
end;

% Stability region,
if nargin < 6,
  H = [];
end;
if ~isempty(H),
 if ~isa(H,'double'),
   error('Invalid stability matrix');
 elseif ~all(size(H)==[2 2]),
   error('Invalid stability matrix');
 end;
 H = (H+H')/2;
else,
  % stability region
 switch symb,
  case {'s','p'}
   H = [0 1;1 0]; 
  case {'z^-1','d'}
   H = [1 0;0 -1];
  case {'z','q'}
   H = [-1 0;0 1];
  otherwise
   error('Invalid input polynomials');
 end;
end;

% Check stability of characteristic polynomial

rd = roots(D);
for i = 1:length(rd),
  if H(1,1)+H(2,1)*rd(i)+H(1,2)*rd(i)'+H(2,2)*rd(i)*rd(i)' >= 0,
    error('Unstable central polynomial');
  end;
end;
  
% PID ?

wpid = nargin > 6;
  
% Margin for ensuring SPRness

mingamma = tol;

% Controller order

xdeg = min(d-deg(a0),d-deg(b0));
ydeg = xdeg;

if xdeg < 0, error(['Degree of central polynomial less than degree of' ...
		 ' plant.']); end;
  
% Build LMI

Na = {}; Nb = {};
na = [a0*eye(n-1) -a1; -a1.' a0];
na = na{:};
for i=1:xdeg+1,
  Na{i} = [zeros(n,(i-1)*n), na, zeros(n,(xdeg-(i-1))*n)];
end

nb = [b0*eye(n-1) -b1; -b1.' b0];
nb = [nb{:} zeros(n,(deg(a0)-deg(b0))*n)];

for i=1:ydeg+1,
  Nb{i} = [zeros(n,(i-1)*n), nb, zeros(n,(xdeg-(i-1))*n)];
end

D = kron(D{:}, eye(n));

if verbose,
  disp('Build LMI');
end;

lmi = sdmpb;

% projection matrix
PI = [eye(d*n) zeros(d*n,n); zeros(d*n,n) eye(d*n)];

[lmi,varQ] = sdmvar(lmi,d*n,'s','Q');
for i=1:xdeg+1,
  [lmi,varx{i}] = sdmvar(lmi,1,1,'x_i');
end
for i=1:ydeg+1,
  [lmi,vary{i}] = sdmvar(lmi,1,1,'y_i');
end
[lmi,vargamma] = sdmvar(lmi,1,1,'gamma');
[lmi,varxy] = sdmvar(lmi,[varx{:} vary{:}],'st','controller parameters');
[lmi,varr] = sdmvar(lmi,1,1,'minnorm');

[lmi,index] = sdmlmi(lmi,(d+1)*n,'H(Q)');
lmi = sdmineq(lmi,index,varQ,PI',[],H);
lmi = sdmineq(lmi,index,vargamma,D');

% determine whether the fixed denominator is demanded (PI or PID controller)
if wpid==1,
  lmi = sdmineq(lmi,-index,0,D',Na{2});
else
  for i=1:xdeg+1,
    lmi = sdmineq(lmi,-index,varx{i},D',Na{i});
  end
end

for i=1:ydeg+1,
  lmi = sdmineq(lmi,-index,vary{i},D',Nb{i});
end

% to avoid the zero solution, we set the coefficient by the highest
%  power of the denominator to 1
% in the case of PI or PID controller, use trace(Q)=1 constraint
if wpid==0,
  [lmi, index3] = sdmlmi (lmi, 1, 'x{n} <= 1');
  lmi = sdmineq(lmi, index3, varx{xdeg+1}, 1);
  lmi = sdmineq(lmi, -index3, 0, 1, .5);
  [lmi, index4] = sdmlmi (lmi, 1, 'x{n} >= 1');
  lmi = sdmineq(lmi, -index4, varx{xdeg+1}, 1);
  lmi = sdmineq(lmi, index4, 0, 1, .5);
else
  [lmi, index_trQ] = sdmlmi(lmi, 1, 'tr Q <= 1');
  for i=1:d*n,
    id=[zeros(1,i-1), 1, zeros(1, d*n-i)];
    lmi = sdmineq(lmi, index_trQ, varQ, id);
  end
  lmi = sdmineq(lmi, -index_trQ, 0, 1, 0.5);
  [lmi, index_trQ2] = sdmlmi(lmi, 1, 'tr Q >= 1');
  for i=1:d*n,
    id=[zeros(1,i-1), 1, zeros(1, d*n-i)];
    lmi = sdmineq(lmi, -index_trQ2, varQ, id);
  end
  lmi = sdmineq(lmi, index_trQ2, 0, 1, 0.5);
end

% minimise the 2-norm of the controller parameters vector
% min r s.t. x'*x <= r^2 (minimum 2-norm) is equivalent to the LMI
% min r s.t. [eye(length(x)) x;x' r^2]>=0

varxy_l = xdeg+ydeg+2;
[lmi, index_minnorm] = sdmlmi(lmi, xdeg+ydeg+3, 'minnorm');
lmi = sdmineq(lmi, -index_minnorm, 0, [eye(varxy_l); zeros(1, varxy_l)], ...
              [eye(varxy_l), zeros(varxy_l,1)]);
lmi = sdmineq(lmi, -index_minnorm, varr, [zeros(xdeg+ydeg+2,1);1]);
for i = 1:xdeg+1,
  L = [zeros(i-1,1); 1; zeros(varxy_l - i + 1, 1)];
  R = [zeros(1,varxy_l), 1];
  lmi = sdmineq(lmi, -index_minnorm, varx{i}, L, R);
end
for i = 1:ydeg+1,
  L = [zeros(xdeg+1+i-1,1); 1; zeros(varxy_l - i + 1 - (xdeg+1), 1)];
  R = [zeros(1,varxy_l), 1];
  lmi = sdmineq(lmi, -index_minnorm, vary{i}, L, R);
end

[lmi, index_gammagrzero] = sdmlmi(lmi, 1, 'gamma > 0');
lmi = sdmineq(lmi, -index_gammagrzero, vargamma, 1, .5);
lmi = sdmineq(lmi, index_gammagrzero, 0, 1, mingamma/2);

% minimize the 2-norm of the controller parameters
lmi = sdmobj(lmi, varr, -1, 1, 'min r');

if verbose,
  disp('Minimize norm of controller coefficients');
  disp('Solve LMI');
end;

pars.fid = verbose;
lmi = sdmsol(lmi,pars);
answer = get(lmi,'feas');

% Retrieve controller

if (answer == 0) | (answer == 1),

  if verbose,
    disp('Problem is feasible');
  end;
 
  Q = get(lmi,'varvalue',varQ);
  for i=1:xdeg+1,
    x{i}=get(lmi,'varvalue',varx{i});
  end
  for i=1:ydeg+1,
    y{i} = get(lmi,'varvalue',vary{i});
  end
  gamma = get(lmi,'varvalue',vargamma);
  r = get(lmi,'varvalue',varr);

  % make a polynomial out of it's coefficients
  xc = zeros(1,xdeg+1);
  if wpid == 1,
    xc(2) = 1;
  else
    for i = 0:xdeg,
      xc(i+1) = x{i+1};
    end;
  end;
  yc = zeros(1,ydeg+1);
  for i = 0:ydeg,
    yc(i+1) = y{i+1};
  end;
  x = pol(xc,xdeg,symb);
  y = pol(yc,ydeg,symb);
 
else

  disp('Problem is infeasible');

  x = []; y = [];
  
end;

% just some testing stuff
%if any(eig(Q)<=0),
%  disp('! Warning ! Q < 0');
%end
%if (trace(Q) >= 1.0001 | trace(Q) <= 0.9999) & wpid==1,
%  disp('trace(Q) != 1');
%end
%if (x{xdeg+1} >= 1.0001 | x{xdeg+1} <= 0.9999) & wpid==0,
%  disp(' x(n) != 1'); 
%end

