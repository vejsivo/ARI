function G = gram(Q,R, arg1, arg2);
%GRAM  Gramian of a stable polynomial matrix fraction
%
% The commands
%    GRAM(N,D)
%    GRAM(N,D,'l') 
% compute the Gramian of the stable rational matrix D\N with 
% N and D left coprime.
% The command
%    GRAM(N,D,'r') 
% computes the Gramian of the stable rational matrix N/D with 
% N and D right coprime.
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.
%
% The controllability Gramian of the state space system (A,B,C,D)
% may be computed as GRAM(B, VAR*I - A, 'l'), while its observability 
% Gramian equals GRAM(C, VAR*I - A, 'r').  VAR is the variable,
% and I stands for the identity matrix with the same size as A.

%       Author(s): M. Hromcik, M. Sebek 22-10-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 28-Oct-1998 10:28:34   $

global PGLOBAL;
opt = 'l';
tol = PGLOBAL.ZEROING;

switch nargin,
  case 1,
    error('Not enough input arguments.');
  case 2,
  case 3,
    if isnumeric(arg1),
      tol = arg1;
    elseif ischar(arg1),
      opt = arg1;
    else error('Invalid command option.');
    end;  
  case 4,
    if isnumeric(arg1),
      tol = arg1;
    elseif ischar(arg1),
      opt = arg1;
    else error('Invalid command option.');
    end;  
    
    if isnumeric(arg2),
      tol = arg2;
    elseif ischar(arg2),
      opt = arg2;
    else error('Invalid command option.');
    end;  
 
end;  	% switch ...

if ~all(size(tol)==1) | ~isreal(tol) | tol<0 | tol>1,
   error('Invalid tolerance.');
end;

eval('Q = pol(Q);', 'error(peel(lasterr));');
eval('R = pol(R);', 'error(peel(lasterr));');
Qv = Q.v; Rv = R.v;

% Var:
if ~strcmp(Qv,Rv),
  if isempty(Qv) | isempty(Rv),
    var = [Qv,Rv];
  else
    warning('Inconsistent variables.');
    var = PGLOBAL.VARIABLE;
  end;
else
  var = Qv;
end; 
Q.v = var; R.v = var;

% Sampling period:
[th,h,Q,R] = testhp(Q,R,var);
if ~th, warning('Inconsistent sampling periods.');
end;

% Denominator:
if size(R,1)~=size(R,2),
   error('Denominator matrix is not square.');
end;

switch opt,
  case 'r',
    % Size:
    if size(Q,2) ~= size(R,2),
       error('Matrices of inconsistent dimensions.');
    end;
    if size(R,2)==0,
       G = []; return;
    end;
  
    % Q,R must be right coprime:
    if ~isprime([Q;R],'r'), error('Matrices are not right coprime.'); end;
    
    % R must be col red
    if rank(lcoef(R,'col')) < min(R.s),
      [R,aux,U] = colred(R);
      Q = Q*U;
    end 
    if any(strcmp(var,{'s','p'})),
      G = sub_r_option_cont(Q,R,tol);
    else
      G = sub_r_option_disc(Q,R,tol,var,opt);
    end;    
  case 'l',
    %[Q,R] = lmf2rmf(Q,R,tol)
    %G = sub_r_option(Q,R,tol);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    	% Size:
    	if size(Q,1) ~= size(R,1),
    	  error('Matrics of inconsistent dimensions.');
      end;
      if size(R,1)==0,
         G = []; return;
      end;
        
      if ~isprime([Q R],'l'),
        error('Matrices are not left coprime.');
      end;

    	% R must be row red
    	if rank(lcoef(R,'row')) < min(R.s),
    	  [R,aux,U] = rowred(R);
    	  Q = U*Q;
    	end  
    	if any(strcmp(var,{'s','p'})),
    	  G = sub_r_option_cont(Q.',R.',tol);
    	else
    	  G = sub_r_option_disc(Q.',R.',tol,var,opt);
    	end;    
    	G = G.';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  otherwise error('Invalid command option.');
end;  	% switch ...

% *************************************************** %
function G = sub_r_option_cont (Q,R,tol);
% SUB_R_OPTION_CONT  Sub-function for 'r' option, continuous.
%
%   SUB_R_OPT(Q,R,tol) computes the gramian of T = Q/R.
%

  % Test if stable:
  if ~ishurwitz(R),
    error('Denominator matrix is not stable.');
  end;

% STEP 1:
degcol = deg(R,'col')-1;
degcol(degcol<0) = 0;
C = axxab(R,Q'*Q,degcol);
if isnan(C), error('The computation of the Gramian failed.'); end;
%if ~isproper(C,R,'r','strictly'), 
%  [aux,C] = rdiv(C,R);
%end; 

% STEP 2:

Rl = lcoef(R,'col');
% Cl:
Cs = C.s;
Cl = zeros(Cs);
c = deg(R, 'col');
for i=1:size(R,2),
  %Cl(:,i) = C{c(i)-1}(:,i);
  Cc = C.c; 
  %eval('Cl(:,i) = C{c(i)-1}(:,i);', 'Cl(:,i) = zeros(Cs(1),1);');
  %eval('Cl(:,i) = C(:,i,c(i));', 'Cl(:,i) = zeros(Cs(1),1);');
  Cl(:,i) = Cc(:,i,c(i));
end;

G = Cl/Rl;
if G(1,1)<0, G = -G; end;

% end .. sub_r_option_cont

% ****************************************************** %
function G = sub_r_option_disc(N,D,tol,var,opt);
% SUB_R_OPTION_DISC  Sub-function for 'r' option, discrete
%
%   SUB_R_OPT(Q,R,tol) computes the gramian of T = Q/R.
%

Ds = D.s;

  % Test if stable:
  if ~isschur(D),
    error('Denominator matrix is not stable.');
  end;  
  
  if any(strcmp(var,{'d','z^-1'})),
    [N,D] = reverse(N,D,tol,opt);
    N.v = 'z';
    D.v = 'z';
  end;

  % Make proper:
  while ~isproper(N,D,'r'),
    D = shift(D,1);
  end;  

  if rank(lcoef(D,'col')) < Ds(1),
     [D,aux,U] = colred(D,tol);
     N = N*U;
  end;

  P = axxab(D,N'*N); 
  if ~isproper(P,D,'r'), 
    [aux,C] = rdiv(P,D);
  end; 
      
  Dl = lcoef(D,'col');
  % Pl:
  Ps = P.s;
  Pl = zeros(Ps);
  c = deg(D, 'col');
  for i=1:size(D,2),
    Pc = P.c; 
    %eval('Pl(:,i) = P(:,i,c(i));', 'Pl(:,i) = zeros(Ps(1),1);');
    Pl(:,i) = Pc(:,i,c(i)+1);
  end; 
  G = Pl/Dl;
  if G(1,1)<0, G = -G; end;
  G = G + G';

% end .. sub_r_option_disc

% References: H. Kwakernaak, Polynomial Computation of Hankel Singular Values, 
%	      CDC 92.

%end .. @pol/gram
