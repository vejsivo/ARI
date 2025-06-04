function [Ak,Bk,Ck,Dk,Ek,clpoles,rho] = ...
                  dsshinf(A,B,C,D,E,nmeas,ncon,gamma,p9,p10,p11,p12,p13)
%DSSHINF  Computation of H-infinity suboptimal compensators for a descriptor standard plant
% 
% The commands
%    [Ak,Bk,Ck,Dk,Ek,clpoles,rho] = DSSHINF(A,B,C,D,E,nmeas,ncon,gamma)
%    [Ak,Bk,Ck,Dk,Ek,clpoles,rho] = DSSHINF(A,B,C,D,E,nmeas,ncon,gamma[,options][,tol])
% compute the suboptimal compensator for the standard H-infinity problem 
% with generalized plant G(s) = C(sE-A)^{-1}B + D  with nmeas measured 
% outputs and ncon control inputs according to the level gamma. Without 
% the option 'all' the routine computes the 'central' compensator 
%    K(s) = Ck*(s*Ek-Ak)^{-1}*Bk
% If the option 'all' is included then after partitioning 
%    Bk = [Bk1 Bk2], Ck = [Ck1        Dk = [0  0
%                          Ck2]             0 Dk22]
% where Bk1 has nmeas columns and Ck1 has ncon rows, all suboptimal 
% compensators may be parametrized as
%    Ek*dxk/dt = Ak*xk  + Bk1*y + Bk2*p
%         u    = Ck1*xk         +   p
%         q    = Ck2*xk +   y   + Dk22*p
%         p    =   U*q
% U is any stable LTI system whose transfer matrix has infinity-norm less 
% than or equal to 1. The central compensator follows by setting U = 0.
%
% The array clpoles contains the finite closed-loop poles obtained for the
% central compensator, that is, the closed poles with size less than 1/tol.
% The default value of the optional input parameter tol is 1e-6. 
%
% The output parameter rho is the smallest singular value of a matrix
% z that arises in the algorithm. For type 2 optimal solutions this
% matrix z is singular and, hence, rho = 0. The parameter may be used
% in searching for the optimal solution.
%
% If gamma is too small so that no suboptimal compensator exists then
% the routine exits with all output arguments empty. In verbose mode 
% a warning message is issued.
%
% Conditions on the input data: If D is partitioned as
%   D = [D11 D12
%        D21 D22]
% where D12 has ncon columns and D21 has nmeas rows then D12 needs
% to have full column rank and D21 full row rank. Use the command DSSREG
% to 'regularize' the system if this condition is not met.
%
% The available options (none or several may be included) besides 'all' are:
%    'notest'    the rank tests on D12 and D21 are omitted
%    'nopoles'   the closed loop poles are not computed, the output argument
%                clpoles is empty

%     Author: H. Kwakernaak, 1997-1998
%     Copyright 1998 PolyX Ltd.
%     Modified by J. Jezek, Aug 2001,  arg checking

% Checks

if nargin < 8
   error('Not enough input arguments.');
end;
if ~isnumeric(A) | ndims(A)>2 | ~isnumeric(B) | ndims(B)>2 | ...
   ~isnumeric(C) | ndims(C)>2 | ~isnumeric(D) | ndims(D)>2 | ...
   ~isnumeric(E) | ndims(E)>2
   error('Invalid 1st - 5th argument; must be numerical matrices.');
end;
[nE,mE] = size(E); [nA,mA] = size(A); [nB,mB] = size(B); 
[nC,mC] = size(C); [nD,mD] = size(D);
if nE~=mE | nA~=mA | nE~=nA | mC~=nA | mA~=nB | nD~=nC | mD~=mB
   error('Matrices of inconsistent dimensions.')
end

if ~isnumeric(nmeas) | length(nmeas)~=1 | ~isreal(nmeas) | ...
   ~isnumeric(ncon)  | length(ncon)~=1  | ~isreal(ncon)
   error('Invalid 6th or 7th argument.');
end;
if ~isnumeric(gamma) | length(gamma)~=1 | ~isreal(gamma)
   error('Invalid 8th argument.');
end;

allcomp = 0; notest = 0; nopoles = 0; tol = 1e-6;
if nargin >= 9
   if ischar(p9)
      if strcmp(p9,'all')
         allcomp = 1;
      elseif strcmp(p9,'notest'); 
         notest = 1;
      elseif strcmp(p9,'nopoles');
         nopoles = 1;
      else
         error('Invalid command option.');
      end
   elseif isnumeric(p9)
      tol = p9;
   else
      error('Invalid 9th argument.');
   end
end
if nargin >= 10
   if ischar(p10)
      if strcmp(p10,'all')
         allcomp = 1;
      elseif strcmp(p10,'notest'); 
         notest = 1;
      elseif strcmp(p10,'nopoles');
         nopoles = 1;
      else
         error('Invalid command option.');
      end
   elseif isnumeric(p10)
      tol = p10;
   else
      error('Invalid 10th argument.');
   end
end
if nargin >= 11
   if ischar(p11)
      if strcmp(p11,'all')
         allcomp = 1;
      elseif strcmp(p11,'notest'); 
         notest = 1;
      elseif strcmp(p11,'nopoles');
         nopoles = 1;
      else
         error('Invalid command option.');
      end
   elseif isnumeric(p11)
      tol = p11;
   else
      error('Invalid 11th argument.');
   end
end
if nargin == 12
   if isnumeric(p12)
      tol = p12;
   else
      error('Invalid 12th argument.');
   end
end

if length(tol)~=1 | ~isreal(tol) | tol<0 | tol>1
   error('Invalid tolerance.');
end

% Initialization

global PGLOBAL;
eval('PGLOBAL.FORMAT;', 'painit;');
verbose = strcmp(PGLOBAL.VERBOSE, 'yes');

k2 = ncon; m2 = nmeas; k1 = mB-k2; m1 = nC-m2; p = nE;
B1 = B(:,1:k1); B2 = B(:,k1+1:k1+k2);
C1 = C(1:m1,:); C2 = C(m1+1:m1+m2,:);
D11 = D(1:m1,1:k1); D12 = D(1:m1,k1+1:k1+k2);
D21 = D(m1+1:m1+m2,1:k1); D22 = D(m1+1:m1+m2,k1+1:k1+k2);


% Further checks

if m1<k2
   error('D12 is not tall.');
elseif m2>k1
   error('D21 is not wide.');
elseif ~notest
   if rank(D12)~=ncon
      error('D12 does not have full column rank.')
   elseif rank(D21)~=nmeas
      error('D21 does not have full row rank.')
   end
end

if gamma == 0
   if verbose
      disp('dsshinf: gamma = 0 is too small.')
   end
   Ak = []; Bk = []; Ck = []; Dk = []; Ek = []; clpoles = []; rho = [];
   return
end
M = D11'*D11/gamma^2; M = eye(size(M))-M;
if rank(M) < length(M)
   if verbose
      disp('dsshinf: gamma is too small.')
   end
   Ak = []; Bk = []; Ck = []; Dk = []; Ek = []; clpoles = []; rho = [];
   return
end


% Form the matrices needed to define the denominator and
% numerator pencils

e = E;
a = A +B1/(eye(k1)-D11'*D11/gamma^2)*D11'*C1/gamma^2;
b2 = B2+B1/(eye(k1)-D11'*D11/gamma^2)*D11'*D12/gamma^2;
c2 = C2+D21/(eye(k1)-D11'*D11/gamma^2)*D11'*C1/gamma^2;
f12 = C1'/(eye(m1)-D11*D11'/gamma^2)*D12;
f21 = D21/(eye(k1)-D11'*D11/gamma^2)*B1';
d22 = D22+D21/(eye(k1)-D11'*D11/gamma^2)*D11'*D12/gamma^2;
r0 = D12'/(eye(m1)-D11*D11'/gamma^2)*D12;
r1 = B1/(eye(k1)-D11'*D11/gamma^2)*B1';
r2 = C1'/(eye(m1)-D11*D11'/gamma^2)*C1;
r3 = D21/(eye(k1)-D11'*D11/gamma^2)*D21';


% Define the denominator and numerator pencils

DenPen = [    -r1+f21'/r3*f21     s*e-a+f21'/r3*c2
           -s*e'-a'+c2'/r3*f21  -r2/gamma^2+c2'/r3*c2 ];
NumPen = [  -r1/gamma^2+b2/r0*b2'  s*e-a+b2/r0*f12'
   -s*e'-a'+f12/r0*b2'    -r2+f12/r0*f12'  ];
DenPen = 0.5*(DenPen+DenPen');
NumPen = 0.5*(NumPen+NumPen');


% Compute the Clements forms of the two pencils

[CleDen,v,fl1] = clements((DenPen.')',0); CleDen = (CleDen.')'; 
[CleNum,w,fl2] = clements((NumPen),0); 
if fl1 == -Inf | fl2 == -Inf % roots on the imaginary axis
   if verbose
      disp('dsshinf: gamma is too small');
   end
   Ak = []; Bk = []; Ck = []; Dk = []; Ek = []; clpoles = []; rho = [];
   return
end
CleDen13 = CleDen(1:p,p+1:2*p);


% Solutions of the GAREs

v11 = v(1:p,1:p);     v12 = v(1:p,p+1:2*p);
v21 = v(p+1:2*p,1:p); v22 = v(p+1:2*p,p+1:2*p);
q1 = v11*f21'+v12*c2';
q2 = v11*b2+v12*f12/gamma^2;
w11 = w(1:p,1:p);     w12 = w(1:p,p+1:2*p);
w21 = w(p+1:2*p,1:p); w22 = w(p+1:2*p,p+1:2*p);
p1 = f12'*w12'+b2'*w11';
p2 = c2*w12'+f21*w11'/gamma^2;


% Compute the compensator generator

epsilon = CleDen13{1};
alpha = -CleDen13{0};
z = v22*w12'+v21*w11'/gamma^2;
Bk1 = q1/r3;
Bk2 = q2-q1/r3*d22;
Ck1 = -r0\p1;
Ck2 = -p2-d22/r0*p1;
Ek = epsilon*z;
Ak = alpha*z-(q2-q1/r3*d22)/r0*p1;


% Compute the closed-loop poles

if nopoles
   clpoles = [];
else
   Ec = [     E      zeros(p,p)
          zeros(p,p)     Ek     ];
   Ac = [    A       B2*Ck1
          Bk1*C2 Ak+Bk1*D22*Ck1 ];
   clpoles = eig(Ac,Ec); 
   clpoles = clpoles(find(~isnan(abs(clpoles))));
   clpoles = clpoles(find(abs(clpoles)<1/tol));
end


% Complete the output

rho = min(svd(z));
if allcomp
   Bk = [Bk1 Bk2];
   Ck = [Ck1; Ck2];
   Dk = [zeros(k2,m2) zeros(k2,k2)
         zeros(m2,k2)    d22      ];
else
   Bk = Bk1;
   Ck = Ck1;
   Dk = zeros(k2,m2);
end

%end .. dsshinf
