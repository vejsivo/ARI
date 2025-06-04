function [Ak,Bk,Ck,Dk,Ek] = dssh2(A,B,C,D,E,nmeas,ncon,tol)
%DSSH2 Descriptor solution of the H2 problem
%
% The command
%    [Ak,Bk,Ck,Dk,Ek] = DSSH2(A,B,C,D,E,nmeas,ncon,tol)
% solves the H2 optimization problem for the standard plant
%    G(s) = C(sE-A)^{-1}B + D  
% with nmeas measured outputs and ncon control inputs.
% The optimal compensator is given by
%    K(s) = Ck*(s*Ek-Ak)^{-1}*Bk+Dk
% The optional parameter tol is a tolerance with default 
% value 1e-10.
%
% Conditions on the input data: If D is partitioned as
%   D = [D11 D12
%        D21 D22]
% where D12 has ncon columns and D21 has nmeas rows then D12 needs
% to have full column rank and D21 full row rank. D22 should be
% the zero matrix. Use the command DSSREG with the option 'D22'
% to 'regularize' the system if these conditions are not met.

%     Author: H. Kwakernaak, April-June 2000
%     Copyright 2000 by PolyX Ltd.
%     Modified by J. Jezek, Aug 2001,  arg checking


% Initialization

if nargin < 7
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

if nargin < 8
   tol = 1e-10;
else
   if ~isnumeric(tol) | length(tol)~=1 | ~isreal(tol) | ...
         tol<0 | tol>1,
      error('Invalid tolerance.');
   end;
end


% Preparation

n = length(A);
k = size(B,2); k2 = ncon; k1 = k-k2;
m = size(C,1); m2 = nmeas; m1 = m-m2;
if n == 0 | k == 0 | m == 0 | k1<=0 | k2<=0 | m1<=0 | m2<=0
   error('Inconsistent dimensions; zero, negative or too high.')
end

B1 = B(:,1:k1);
B2 = B(:,k1+1:k);
C1 = C(1:m1,:);
C2 = C(m1+1:m,:);
D11 = D(1:m1,1:k1);  D12 = D(1:m1,k1+1:k);
D21 = D(m1+1:m,1:k1); D22 = D(m1+1:m,k1+1:k);
Q = D12'*D12;
R = D21*D21';


% Checks

if issingular(Q) | issingular(R)
   disp('D12 or D21 does not have full rank.')
   error('Use dssreg to regularize the plant.')
elseif norm(D22,1)>tol*norm([A B;C D],1)
   disp('D22 is nonzero.')
   error('Use dssreg with option D22.')
end
   

% Solve the regulator and observer GAREs 

[X,F] = gare(A,B2,C1,D12,E,eye(m1),Q);
[Y,K] = gare(A',C2',B1',D21',E',eye(k1),R); Y = Y'; K = K';


% Check the roots

rts1 = eig(A-B2*F,E); rts1 = rts1(find(abs(rts1)<1/sqrt(tol)));
rts2 = eig(A-K*C2,E); rts2 = rts2(find(abs(rts2)<1/sqrt(tol)));
if any(real(rts1)>sqrt(tol))
   error('The plant has unstable fixed poles.')
elseif any(real(rts2)>sqrt(tol))
   error('The plant has unstable fixed poles.')
end


% Initiate the final computations

[sEA2,rk,U] = colred(pzer(s*E-A+B2*F));
sEA2 = pzer(sEA2); CC2 = pzer((C1-D12*F)*U);
[n1,d1] = rmf2lmf(CC2,sEA2);
[sEA1,rk,U] = rowred(pzer(s*E-A+K*C2));
sEA1 = pzer(sEA1); BB1 = pzer(U*(B1-K*D21));
[n2,d2] = lmf2rmf(BB1,sEA1);
nb12 = d1*D12+n1*B2; db12 = d1;
nb21 = C2*n2+D21*d2; db21 = d2;
[n12,d12] = lmf2rmf(nb12,db12);
[n21,d21] = rmf2lmf(nb21,db21);


% Existence checks

if k2<m1
   k12 = null(n12');
   [nn1,dd1] = r2l(k12',d1);
   if any(abs(real(roots(dd1)))<tol)
      error('The plant has marginally stable fixed poles that cannot be cancelled.')
   end
   pp1 = ldiv(nn1*n1*B1,dd1)+k12'*D11;
   if any(deg(pp1,'row')>=deg(k12','row'))
      error('The closed-loop transfer matrix cannot be made strictly proper.')
   end
elseif m2<k1
   k21 = null(n21)';
   [nn2,dd2] = l2r(k21',d2);
   if any(abs(real(roots(dd2)))<tol)
      error('The plant has marginally stable fixed poles that cannot be cancelled.')
   end
   pp2 = rdiv(C1*n2*nn2,dd2)+D11*k21';
   if any(deg(pp2,'col')>=deg(k21','col'))
      error('The closed-loop transfer matrix cannot be made strictly proper.')
   end
end


% Prepare the computation of L

[mb21,dbo] = rmf2lmf(n21,d2');
[m12,do] = lmf2rmf(n12,d1');
p = n1*B1*d2+n1*B2*F*n2+d1*D12*F*n2+d1*D11*d2;
pb = mb21*p'*m12;


% Solve the Diophantine equations

[c2,c1] = axybc(dbo*d21,d12*do,pb);
[a2,a1] = axybc(d21,dbo,eye(size(dbo)));
[b2,b1] = axybc(do,d12,eye(size(do)));


% Compute L

P = ldiv(a1*c1,d21)+rdiv(c2*b2,d12);
[nl,dl] = l2r([c1; b1],[dbo zeros(m2,k2); zeros(k2,m2) do]);
nl = P*dl+[a2 c2]*nl;
nl = R\nl;
dl = Q*dl;
[nl,dl] = r2l(nl,dl);
[alpha,beta,gamma,delta,epsilon] = rmf2dss(nl',dl');


% Compute the compensator

nn = length(alpha);
Ek = [E zeros(n,nn); zeros(nn,n) epsilon];
Ak = [A-B2*F-K*C2+B2*delta*C2 -B2*gamma; -beta*C2 alpha];
Bk = [K-B2*delta; beta];
Ck = [-F+delta*C2 -gamma];
Dk = -delta;


% Auxiliary routines

function [NR,DR] = l2r(NL,DL)
% Left to right conversion with precautions for nonproper fractions
[DL,rk,U] = rowred(DL); NL = U*NL;
[Q,R] = ldiv(NL,DL);
[R,DR] = lmf2rmf(R,DL);
NR = pzer(Q*DR+R);
return

function [NL,DL] = r2l(NR,DR)
% Right to left conversion with precautions for nonproper fractions
[NL,DL] = l2r(NR',DR');
NL = NL'; DL = DL';
return

%end .. dssh2
