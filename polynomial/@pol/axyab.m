function [x,y] = axyab(a,blp,brp,offset,varargin)
%AXYAB   Non-symmetric polynomial equation solver
%
% For discrete-time polynomials A,B, the command
%    [X,Y] = AXYAB(A,BL,BR,OFFSET) 
% solves the scalar polynomial equation
%    A'X + Y'A = BL' + SHIFT(BR, OFFSET).
% A is a stable polynomial and it is assumed that 
%    MAX(DEG(BL), DEG(BR)) <= DEG(A).
% The solution polynomials X, Y are such that 
%    MAX(DEG(X), DEG(Y)) <= DEG(A).
%
% For continuous-time polynomials A,B, the argument OFFSET
% plays no role and, if it is the last argument, it may be
% omitted. It is assumed OFFSET=0, i.e. the equation to be
% soled is    A'X + Y'A = BL' + BR .
%
% See also TSP/AXYAB. That routine solves the scalar polyomial
% equation  A'X + Y'A = B  with a scalar polynomial A in variable
% 'z' or 'z^-1' and a two-sided polynomial B.
%
% The command
%    AXYAB(A,BL,BR,OFFSET,'syl') 
% uses the Sylvester matrix algorithm. This is the default method.
%
% The command
%    AXYAB(A,BL,BR,OFFSET,'red') 
% uses the polynomial reduction algorithm, a version of the Euclidean 
% division algorithm for polynomials.
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.
%
% See also: POL/AXXAB, POLPART.

%    Author: D. Henrion, January 19, 1999.
%    Copyright 1998-99 by Polyx, Ltd.

global PGLOBAL;

eval('PGLOBAL.FORMAT;','painit;');

if nargin < 3,
 error('Not enough input arguments.');
end;

eval('a = pol(a); blp = pol(blp); brp = pol(brp);', ...
   'error(peel(lasterr));');

[tv,typec,a,blp,brp] = testvp3cd(a,blp,brp);
if tv==0,
   warning('Inconsistent variables.');
elseif tv==-1,
   error('Inconsistent variables, continuous-time versus discrete-time.');
end;

[th,Xh,a,blp,brp] = testhp3(a,blp,brp,typec);
if th==0,
   warning('Inconsistent sampling periods.');
end;

if strcmp(typec, 's') | strcmp(typec, 'p'),
   type = 'continuous';
else
   type = 'discrete';
   if nargin<4,
      error('Not enough input arguments.');
   end;
   if ~isa(offset,'double') | length(offset)~=1 | ...
      ~isreal(offset) | floor(offset)~=offset,
         error('Invalid offset.');
   end;
   offset = max(offset,0);
end;

if strcmp(type, 'continuous')
   b = blp'+brp;
else
   b = rev(conj(blp.'))+shift(brp,offset);
end;

if any(any(isnan(a))) | any(any(isnan(b))),
 error('Polynomial is not finite.');
elseif any(any(isinf(a))) | any(any(isinf(b))),
 error('Polynomial is not finite.');
end;

if ~isreal(a) | ~isreal(b),
 error('Complex polynomials not yet implemented.');
end;

% verbose level
verbose = strcmp(PGLOBAL.VERBOSE, 'yes');

[ra, ca] = size(a); dega = a.d;
[rb, cb] = size(b); degb = b.d;

if ra ~= 1 | ca ~= 1 | rb ~= 1 | cb ~= 1,
   error('Invalid arguments; must be scalar polynomials.');
end;

n = ra;

% Default options.

tol = [];
method = [];

% Parse optional input arguments.

invalid = 0;
lv = length(varargin);
if lv>0,
 for i = 1:lv,
  arg = varargin{i};
  if ~isempty(arg),
   if isa(arg, 'char'), % valid strings
    if strcmp(arg, 'syl'),
     method = 'syl';
    elseif strcmp(arg, 'red'),
     method = 'red';
    else,
     error('Invalid command option.');
    end;     
   elseif isa(arg, 'double'), % matrix or scalar
    if ~any(size(arg) - [1 1]), % scalar = tolerance
     tol = arg;
    else,
     error(['Invalid ',nth(i+4),' argument.']);
    end;
   else
    error(['Invalid ',nth(i+4),' argument.']);
   end;
  end;
 end;
end;

% tolerance for zeroing:
% relative to optional input parameter TOL and to elements in A and B
ac = a.c; bc = b.c;
me = max(norm(ac(:,:), 'inf'), norm(bc(:,:), 'inf'));
if isempty(tol),
 tolzero = PGLOBAL.ZEROING * me;
else
 tolzero = tol * me;
end;

% tolerance for rank computation:
tolrank = tolzero;

solution = 1;

if degb < 0,

 % b is zero, so are x and y
 x = 0; degx = 0; 
 y = 0; degy = 0;

elseif strcmp(type, 'continuous'),
 
 % CONTINUOUS-TIME

 if isempty(method), method = 'syl';
 end;

 if strcmp(method, 'red'),

  % ----------------------------------------
  % CONTINUOUS-TIME
  % POLYNOMIAL REDUCTION ALGORITHM
  % ----------------------------------------

  a2 = zeros(1, dega+1); a2(:) = a.c;
  b2 = zeros(1, degb+1); b2(:) = b.c/2;

  % separation of odd and even parts of a and b

  na = floor(dega/2)+1; nb = floor(degb/2)+1; m = max(na, nb);
  ae = zeros(1,m); ao = ae; be = ae; bo = ae;
  a2(1, dega+2) = 0; b2(1, degb+2) = 0;

  for i = 1:na,
    ae(i) = a2(1, (i-1)*2+1);
    ao(i) = a2(1, 2*i);
  end;
  for i = 1:nb,
    be(i) = b2(1, (i-1)*2+1);
    bo(i) = -b2(1, 2*i);
  end;

  % forward construction of coefficients

  coefs = []; nit = 1;

  while length(find(ao)) > 0, % while polynomial ao(s) is not zero
    a1 = ao(1); a2 = ae(1);
    if (abs(a1) < tolzero) | (abs(a2) < tolzero),
     error('Invalid 1st polynomial; not strictly stable.');
    end;
    coefs = [coefs; a2/a1 be(1)/a2 bo(1)/a1];
    be = [be(2:m) - coefs(nit,2)*ae(2:m) 0];
    bo = [bo(2:m) - coefs(nit,3)*ao(2:m) 0];
    aen = ao; ao = [ae(2:m) - coefs(nit,1)*ao(2:m) 0]; ae = aen;
    nit = nit + 1;
  end;

  % backward substitutions

  ue = zeros(1,m); uo = ue; ve = ue; vo = ue;
  ue(1:length(be)) = be/ae(1);
  vo(1:length(bo)) = bo/ae(1);

  for i = nit-1:-1:1,
    uon = ue - coefs(i,1)*uo;
    ue = [coefs(i,2) uo(1:m-1)]; uo = uon;
    ven = [coefs(i,3) vo(1:m-1)] - coefs(i,1)*ve;
    vo = ve; ve = ven;
  end;

  % construction of solutions x(s) and y(s)

  xe = ue + ve; xo = -(uo + vo); x = [];
  ye = ue - ve; yo = -(uo - vo); y = [];

  for i = 1:m,
   x = [x xe(i) xo(i)];
   y = [y ye(i) yo(i)];
  end;

  degx = 2*m-1; degy = 2*m-1;

 elseif strcmp(method, 'syl') % method

  % ---------------------------------
  % CONTINUOUS-TIME
  % SYLVESTER MATRIX METHOD
  % ---------------------------------

  b2 = zeros(1, degb+1); b2(:) = b.c / 2;

  % maximum possible degree for x and y
  degx = dega; degy = dega;
  b2 = [b2 zeros(1,dega+degx-degb)];
  R = [sylv(a', degx); sylv(a, degx)];

  if rank(R, tolrank) ~= size(R, 2),
    xy = 2 * b2 * pinv(R);
    % check
    if norm(xy*R-2*b2) > tolzero,
     solution = 0;
    end;
  else
    xy = 2 * b2 / R;
  end;

  if solution,

   x = xy(1:degx+1); y = xy(degx+2:2*(degx+1));
 
   for i = 2:2:degx+1,
    y(i) = -y(i);
   end;

  end;

 else

  error('Invalid input argument.');

 end; % if method

else % continuous or discrete

 % DISCRETE-TIME

 if isempty(method), method = 'red'; end;

 if strcmp(method, 'syl'),
   warning('AXYAB: Sylvester matrix method not yet implemented in discrete-time.');
   disp('Switch to polyomial reduction method.');
   method = 'red';
 end;

 if strcmp(method, 'red'),

  % ----------------------------
  % DISCRETE-TIME
  % POLYNOMIAL REDUCTION METHOD
  % ----------------------------

  dshift = ceil(degb/2);
  degbl = max(dshift, degb-dshift);

  a2 = zeros(1,dega+1); a2(:) = a.c;
  b2 = zeros(1,degb+1); b2(:) = b.c;

  br = b2(1,dshift+1)/2; bl = br;
  if dshift < degb, br = [br b2(1,dshift+2:degb+1)]; end;
  br = [br zeros(1,degbl-degb+dshift)];
  if dshift > 0, bl = [bl b2(1,dshift:-1:1)]; end;
  bl = [bl zeros(1,degbl-dshift)];

  degx = max([dega degbl]); degy = degx;
  a2 = [a2 zeros(1,degx-dega)]; olda = a2;
  degbr = degbl;

  % forward construction of coefficients
 
  backward = []; nit = 0;
 
  while (dega > 0) | (degbr > 0) | (degbl > 0),
   nit = nit + 1;
   if abs(a2(1)) < tolzero,
    error('The first input polynomial must be strictly stable.');
   end;
   if degbr > dega,
     coef = br(degbr+1)/a2(1);
     for i = 1:degbr,
       br(i) = br(i) - coef*a2(degbr+2-i);
    end;
      br(degbr+1) = 0;
     backward = [backward; 1 coef degbr];
     degbr = degbr - 1;
   elseif degbl > dega,
     coef = bl(degbl+1)/a2(1);
     for i = 1:degbl,
       bl(i) = bl(i) - coef*a2(degbl+2-i);
     end;
     bl(degbl+1) = 0;
     backward = [backward; 2 coef degbl];
     degbl = degbl - 1;
   else % degbl <= dega  and  degbr <= dega
     coef = a2(dega+1)/a2(1);
     for i = 1:dega,
       an(i) = a2(i) - coef*a2(dega+2-i);
     end;
     a2 = [an(1:dega) 0];
     backward = [backward; 3 coef dega];
     dega = dega - 1;
   end; 
 
  end %  while

  % backward substitutions

  x = zeros(1,degx+1); y = x; xn = x; yn = y;
  x(1) = br(1)/a2(1);
  y(1) = bl(1)/a2(1);
  for i = nit:-1:1,
   type = backward(i,1);
   coef = backward(i,2);
   deg = backward(i,3);
   if type == 1,
    x(deg+1) = x(deg+1) + coef;
   elseif type == 2,
     y(deg+1) = y(deg+1) + coef;
   else % type == 3
     for j = 1:deg+1,
       xn(j) = x(j) - coef*y(deg+2-j);
       yn(j) = y(j) - coef*x(deg+2-j);
      end;
     x = xn; y = yn;
   end; % if type
  end; % for i 

  % make constant term of x real
  if abs(imag(x(1))) > tolzero,
   coef = sqrt(-1)*imag(x(1))/real(olda(1));
   for i = 1:degx+1,
    x(i) = x(i) - coef * olda(i);
   end;
  end;

 else

  error('Invalid input argument.');

 end;

end;

if solution,
  x = pzer(pol(x, degx)); pprop(x, typec);
  y = pzer(pol(y, degy)); pprop(y, typec);

  if verbose,
   if strcmp(type, 'continuous'),
    residue = norm(a'*x+y'*a-b);
   else
    acjg = a'; dega = a.d;
    ycjg = y'; degy = y.d;
    bcjg = blp'; degb = offset;
    m = max([dega degy degb]);
    p1 = shift(acjg*x, m - dega); pprop(p1, typec);
    p2 = shift(ycjg*a, m - degy); pprop(p2, typec);
    p3 = shift(bcjg, m - degb); pprop(p3, typec);
    p4 = shift(brp, m); pprop(p4, typec);
    residue = norm(p1+p2-p3-p4);
   end;
   if residue > tolzero,
    warning(['AXYAB: Final check. High residue = ' num2str(residue)]);
   else
    disp(['AXYAB: Final check is OK.']);
   end;
  end;

else

  x = pol(NaN, 0);
  y = pol(NaN, 0);

  if verbose,
   disp('AXYAB: No polynomial solution was found.');
  end;

end;

X.h  = Xh; Y.h = Xh;

%end .. @pol/axyab


