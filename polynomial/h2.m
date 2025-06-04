function [Y,X,clpoles,fixed] = h2(N,D,nmeas,ncon,tol,check)
%H2   H2-optimization
%
% Given a continuous-time generalized plant of the form
%    [ z ]    [ G11  G12 ] [ v ]
%    [   ]  = [          ] [   ]
%    [ y ]    [ G21  G22 ] [ u ]
% with G represented in the left coprime polynomial matrix  
% fraction form G = D\N the command
%    [Y,X,clpoles,fixed] = H2(N,D,nmeas,ncon)
% computes the compensator u = K*y in right polynomial matrix 
% fraction form K = Y/X that minimizes the 2-norm of the 
% closed-loop transfer matrix H from v to z. The input
% parameter ncon is the number of control inputs and nmeas
% is the number of measured outputs.
%
% The output parameter clpoles contains the (non-fixed) 
% closed-loop poles and fixed the fixed plant poles.
%
% In the optional forms
%    [Y,X,clpoles,fixed] = H2(N,D,nmeas,ncon,tol)
%    [Y,X,clpoles,fixed] = H2(N,D,nmeas,ncon,'check')
%    [Y,X,clpoles,fixed] = H2(N,D,nmeas,ncon,tol,'check')
% the parameter tol is an optional tolerance. Its default  
% value is 1e-8.
%
% If the option 'check' is present then the routine checks 
% whether the H2-optimization problem has a solution, and 
% exits if no solution exists. If the option is not invoked 
% then the routine produces a solution even if none exists. 
% In this case the closed-loop transfer matrix either has 
% poles on the imaginary axis or is not strictly proper.

%    Author: H. Kwakernaak, August-September, 1999; May, 2000
%    Copyright 1999, 2000, PolyX Ltd
%    Modified by J. Jezek, August 2001,  arg checking


% Preliminaries

if nargin < 4
   error('Not enough input arguments.');
elseif nargin == 4
   tol = 1e-8;
   check = 'none';
elseif nargin == 5
   if ischar(tol)
      check = tol;
      tol = 1e-8;
   elseif ~isnumeric(tol),
      error('Invalid 5th argument.');
   end
else
   if ~isnumeric(tol)
      error('Invalid 5th argument.');
   end;
   if ~ischar(check)
      error('Invalid 6th argument.');
   end;
end

if length(tol)~=1 | ~isreal(tol) | tol<0 | tol>1,
   error('Invalid tolerance.');
end;

eval('N = pol(N); D = pol(D);', ...
   'error(peel(lasterr));')

[rN, cN] = size(N);
[rD, cD] = size(D);
if rN==0 | cN==0,
   error('Numer of inputs or outputs is zero.');
end;

eval('[N,D] = testdnd(N,D,''l'');', ...
   'error(peel(lasterr));');
Var = ''; h = 0;
eval('[Var,h,N,D] = testvhnd(N,D);', ...
   'error(peel(lasterr));');

if strcmp(Var, 'z^-1') | strcmp(Var, 'd'),
   [N,D] = reverse(N,D,'l',tol);
   Var = 'z'; N.v = 'z'; D.v = 'z';
end;

if ~isa(ncon,'double') | length(ncon)~=1 | ~isreal(ncon) | ...
      ~isa(nmeas,'double') | length(nmeas)~=1 | ~isreal(nmeas) | ...   
      ncon<=0 | nmeas<=0,
   error('Invalid 3rd or 4th argument.');
end;
k2 = ncon; m2 = nmeas;
[m,k] = size(N);
if k2>=k | m2>=m,
   error('Invalid 3rd or 4th argument; too high.');
end;
k1 = k-k2; m1 = m-m2;


% Bring N, D into standard left fractional form

[D1,r,U] = rowred(D(:,1:m1));
%D = [D1 U*D(:,m1+1:m)]; N = U*N;        % J. Jezek:
D = horzcat(D1,U*D(:,m1+1:m)); N = U*N;  %slightly quicker
D11 = D(1:m1,1:m1);     D12 = D(1:m1,m1+1:m);
D21 = D(m1+1:m,1:m1);   D22 = D(m1+1:m,m1+1:m);
N11 = N(1:m1,1:k1);     N12 = N(1:m1,k1+1:k);
N21 = N(m1+1:m,1:k1);   N22 = N(m1+1:m,k1+1:k);
rootsD11 = roots(D11);
[N21,r,U] = rowred(N21);
D22 = U*D22; N22 = U*N22;


% Compute the reduced 22 fractions and divisors

[Nhat22,Dhat22] = l2r(N22,D22);
[Ncheck22,Dcheck22] = r2l(Nhat22,Dhat22);
Do = rdiv(D22,Dcheck22);
rootsDo = roots(Do);
[Nb12,Dbo] = l2r(N12*Dhat22-D12*Nhat22,D11);
[Nb12,r,U] = colred(Nb12); Dbo = Dbo*U;


% Checks

if ~(strcmp(check,'check') | strcmp(check,'none'))
   error('Invalid command option.')
end
%if rank(Nb12)>m1      changed  31-Aug-2001  by J.Jezek
%if rank(Nb12)<m1      changed  15-Sep-2001  by J.Jezek
if rank(Nb12)<k1
   error('G12 does not have full column rank.');
%elseif rank(N21)>k1   changed  31-Aug-2001  by J.Jezek
%elseif rank(N21)<k1   changed  15-Sep-2001  by J.Jezek
elseif rank(N21)<m2
   error('G21 does not have full row rank.')
end
if any(abs(real(roots(N21)))< tol)
   warning('h2: N21 has zeros on the imaginary axis.')
end
if any(abs(real(roots(Nb12)))< tol)
   warning('h2: Nb12 has zeros on the imaginary axis.')
end
if any(real(rootsD11)> tol)
   error('The generalized plant has unstable fixed poles.')
end
if any(real(rootsDo)> tol)
   error('The generalized plant has unstable fixed poles.')
end


% Checks for the existence of a solution (optional)

if strcmp(check,'check')
   
   if m2<k1
      M21 = null(N21); M21 = M21';
      [nr,dr] = l2r(N11*M21',D11);
      if any(real(roots(dr))> -tol)
         error('The plant has fixed poles on the imaginary axis that cannot be canceled.');
      end
      pr = rdiv(nr,dr);
      if ~all(deg(pr,'col') < deg(M21,'row')')
         error('The closed-loop system transfer matrix cannot be made strictly proper.');
      end
   end % first check

   if k2<m1
      [Nb,Db] = l2r(N,D);
      [Db1,r,U] = colred(Db(1:k1,:));
      Db = [Db1; Db(k1+1:k,:)*U]; Nb = Nb*U;
      Db11 = Db(1:k1,1:k1);
      Nb11 = Nb(1:m1,1:k1);
      Nbb12 = colred(Nb(1:m1,k1+1:k));
      Mb12 = null(Nbb12');
      [nl,dl] = r2l(Mb12'*Nb11,Db11);
      if any(real(roots(dl))> -tol)
         error('The plant has fixed poles on the imaginary axis that cannot be canceled.');
      end
      pl = ldiv(nl,dl);
      if ~all(deg(pl,'row')<deg(Mb12,'col')')
         error('The closed-loop system transfer matrix cannot be made strictly proper.');
      end
   end % second check
end % check


% Do the two spectral factorizations:

Phi = spcof(N21*N21'); 
Psi = spf(Nb12'*Nb12);


% Solve the Diophantine equation

[Xo,Yo] = axbyc(Dcheck22,-Ncheck22,eye(m2));


% Solve the two two-sided equations

R1 = [    eye(m1)    zeros(m1,m2)  ];
T =  [     D11      D12*Xo-N12*Yo
         zeros(m2,m1)      Do      ];
R2 = [ N11
       N21 ];        
[B21,B11] = axybc(D11,Phi',N11*N21');
B1 = [B11; Phi];
B2 = [B21; zeros(m2,m2)];
[A21,A11] = axybc(Psi',D11,Nb12');
[A22,A12] = axybc(Psi',Do,-A11*(D12*Xo-N12*Yo));
A1 = [A11 A12];
A2 = [A21 A22];


% Compute the polynomial part

P = ldiv(A1*R2*N21'+Nb12'*R1*B2-A1*T*B2,Psi');
P = rdiv(P,Phi');


% Compute the data for the optimal parameter

[phi,do] = l2r(Phi,Do);
[psi,dbo] = r2l(Psi,Dbo);
b1 = [B11*do; phi];
[D11a,rk,U] = colred(D11);
dboA21 = pzer(dbo*A21*U,tol*norm(dbo*A21,1));
a21 = rdiv(dboA21,D11a);
a2 = [a21 dbo*A22];
q = a2*[eye(m1) -T(1:m1,m1+1:m1+m2); zeros(m2,m1) eye(m2)]*b1;
p = dbo*P*do+q;
[pb,psib] = l2r(p,psi);


% Compute the optimal compensator K = Y/X

r = colred([phi*psib; -pb]);
Y = [Yo Dhat22]*r;
X = [Xo Nhat22]*r;
Y = pzer(Y,tol);
X = pzer(X,tol);
[X,rk,U] = colred(X);
Y = pzer(Y*U,tol);
for i = 1:m2
   nn(i) = max(max(abs(X{:}(:,i))));
   if nn(i) < tol, nn(i) = 1; end
end
Y = Y/diag(nn);
X = X/diag(nn);


% Compute the variable closed-loop poles


clpoles = roots(Dcheck22*X-Ncheck22*Y);
if any(abs(real(clpoles))< tol)
   warning('h2: The closed-loop system has one or more poles on the imaginary axis')
end


% Compute the fixed closed-loop poles

if nargout >= 4
   fixed = [rootsD11; rootsDo];
end


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

%end .. h2
