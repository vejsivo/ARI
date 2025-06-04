function [a,b,c,d,e] = dssreg(A,B,C,D,E,nmeas,ncon,arg8,arg9,arg10)
%DSSREG     "Regularizes" a standard descriptor plant
% 
% The commands
%    [a,b,c,d,e] = DSSREG(A,B,C,D,E,NMEAS,NCON)
%    [a,b,c,d,e] = DSSREG(A,B,C,D,E,NMEAS,NCON[,tol][,option1][,option2])
% transform the generalized plant
%    Edx/dt = Ax + B[w; u]
%    [y; z] = Cx + D[w; u]
% where the dimension of y is NMEAS and the dimension of u is NCON,
% to an equivalent generalized plant 
%    edx/dt = ax + b[w; u]
%    [y; z] = cx + d[w; u]
% with
%    d = [d11 d12
%         d21 d22]
% such that d12 has full column rank and d21 full row rank.
% Equivalence means that the two plants have the same transfer
% matrices.
%
% The optional tolerance parameter tol is used in the various rank 
% tests. It has the default value 1e-12.
%
% Two options may be included. The option 'D11' modifies the 
% representation so that the term with d11 is absent. The option 
% 'D22' removes the term with d22.
%
% In verbose mode the routine displays a relative error based on the 
% differences of the frequency response matrices of the transformed 
% and the original plant at the frequencies 1, 2, ..., 10. In non-
% verbose mode a warning is issued if this error is larger than 1e-6

% Author: H. Kwakernaak, June, 1998. Modified May, 2000.
% Copyright 1998, 2000 by PolyX Ltd.
% Modified by J. Jezek, 13-Aug-2001, arg checking

% Checks

if nargin<7,
   error('Not enough input arguments.');
end;
if ~isnumeric(A) | ~isnumeric(B) | ~isnumeric(C) | ...
   ~isnumeric(D) | ~isnumeric(E) | ndims(A)>2 | ...
   ndims(B)>2 | ndims(C)>2 | ndims(D)>2 | ndims(E)>2,
      error('Invalid 1st - 5th argument; must be numerical matrices.');
end;
[nE,mE] = size(E); [nA,mA] = size(A); [nB,mB] = size(B); 
[nC,mC] = size(C); [nD,mD] = size(D);
if nE~=mE | nA~=mA | nE~=nA | mC~=nA | mA~=nB | nD~=nC | mD~=mB
   error('Matrices of inconsistent dimensions.'); 
end
if ~isa(nmeas,'double') | length(nmeas)~=1 | ~isreal(nmeas) | ...
      nmeas<0 | nmeas>nC,
   error('Invalid 6th argument.');
end;
if ~isa(ncon,'double') | length(ncon)~=1 | ~isreal(ncon) | ...
      ncon<0 | ncon>mB,
   error('Invalid 7th argument.');
end;

% Initialization

global PGLOBAL;
eval('PGLOBAL.FORMAT;', 'painit;');
verbose = strcmp(PGLOBAL.VERBOSE, 'yes');

tol = 1e-12; option1 = ''; option2 = ''; option3 = '';
if nargin >= 8,
   if isa(arg8,'double'), tol = arg8;
   elseif isa(arg8,'char'), option1 = arg8;
   else error('Invalid 8th argument.');
   end;
end;
if nargin >= 9,
   if isa(arg9,'double'), tol = arg9;
   elseif isa(arg9,'char'), option2 = arg9;
   else error('Invalid 9th argument.');
   end;
end;
if nargin == 10,
   if isa(arg10,'double'), tol = arg10;
   elseif isa(arg10,'char'), option3 = arg10;
   else error('Invalid 10th argument.');
   end;
end;

if length(tol)~=1 | ~isreal(tol) | tol<0 | tol>1,
   error('Invalid tolerance.');
end;

k2 = ncon; m2 = nmeas; k1 = mB-k2; m1 = nC-m2; p = nE;
B1 = B(:,1:k1); B2 = B(:,k1+1:k1+k2);
C1 = C(1:m1,:); C2 = C(m1+1:m1+m2,:);
D11 = D(1:m1,1:k1); D12 = D(1:m1,k1+1:k1+k2);
D21 = D(m1+1:m1+m2,1:k1); D22 = D(m1+1:m1+m2,k1+1:k1+k2);

e = E; a = A; b = B; c = C; d = D;
b1 = B1; b2 = B2;
c1 = C1; c2 = C2;
d11 = D11; d12 = D12;
d21 = D21; d22 = D22;


% Rank checks

if m1<k2 
   error('D12 is not tall.'); 
elseif m2>k1
   error('D21 is not wide.');
end
[u1,s1,v1] = svd(D12);
r1 = length(find(diag(s1)>tol*norm([C1 D11 D12],1)));
[u2,s2,v2] = svd(D21);
r2 = length(find(diag(s2)>tol*norm([C2 D21 D22],1)));
[U,S,V] = svd(D11);
r3 = length(find(diag(S)>tol*norm([C1 D11 D12])));
if r1 == k2 & r2 == m2 & r3 == 0, return, end


% Left and right null spaces of E

[N,s,M] = svd(E);
n = length(find(diag(s)>tol*s(1,1)));
N = N(:,n+1:p)';
M = M(:,n+1:p);
q = size(N,1);


% Rank increase of D12 by system transformation

if r1 < m1
   if isempty(M) 
      Y = zeros(q,k2); 
   elseif norm(C1*M) < tol*norm([C1' M])
      Y = zeros(q,k2);
   elseif norm(D12) < tol*norm(C1)
      [U,S,V] = svd(C1*M);
      Y = V*eye(q,k2);
   else
      U1 = u1(:,r1+1:k2);
      Z = null(U1'*C1*M); 
      s1 = size(Z,2);
      if r1+s1 <= k2
         V2 = v1(r1+1:r1+s1,:);
         Y = Z*V2;
      else
         Z1 = Z(:,1:k2-r1);
         V2 = v1(r1+1:k2,:);
         Y = Z1*V2;
      end
   end
   B2 = B2+A*M*Y;
   D12 = D12+C1*M*Y;
   D22 = D22+C2*M*Y;
end


% Rank increase of D21 by system transformation

if r2 < m2
   if isempty(N) 
      X = zeros(m2,q);
   elseif norm(N*B1) < tol*norm([N' B1])
      X = zeros(m2,q);
   elseif norm(D21) < tol*norm(C2)
      [U,S,V] = svd(N*B1);
      X = eye(m2,q)*U';
   else
      V1 = v2(:,r2+1:k1);
      W = null(V1'*B1'*N'); W = W';
      s2 = size(W,1);
      if r2+s2 <= m2
         U2 = u2(:,r2+1:r2+s2);
         X = U2*W;
      else
         W1 = W(1:m2-r2,:);
         U2 = u2(:,r2+1:m2);
         X = U2*W1;
      end
   end
   C2 = C2+X*N*A;
   D21 = D21+X*N*B1;
   D22 = D22+X*N*B2;
end


% Rank completion of D12 by adding pseudo states

[U,S,V] = svd(D12);
r1 = length(find(diag(S)>tol*norm([C1 D11 D12])));
if r1 < k2
   U2 = U(:,r1+1:m1);
   V2 = V(:,r1+1:k2);
   So = eye(m1-r1,k2-r1);
   E = [       E         zeros(p,k2-r1)
        zeros(k2-r1,p) zeros(k2-r1,k2-r1)];
   A = [       A         zeros(p,k2-r1)
        zeros(k2-r1,p)  eye(k2-r1,k2-r1)];
   B1 = [B1; zeros(k2-r1,k1)];
   B2 = [B2; V2'];
   C1 = [C1 U2*So];
   D12 = D12+U2*So*V2';
   C2 = [C2 zeros(m2,k2-r1)];
end


% Rank completion of D21 by adding pseudo states

[U,S,V] = svd(D21);
r2 = length(find(diag(S)>tol*norm([C2 D21 D22])));
p = length(E);
if r2 < m2
   U2 = U(:,r2+1:m2);
   V2 = V(:,r2+1:k1);
   So = eye(m2-r2,k1-r2);
   E = [       E         zeros(p,m2-r2)
        zeros(m2-r2,p) zeros(m2-r2,m2-r2)];
   A = [       A         zeros(p,m2-r2)
        zeros(m2-r2,p)  eye(m2-r2,m2-r2)];
   B1 = [B1; So*V2'];
   B2 = [B2; zeros(m2-r2,k2)];
   C1 = [C1 zeros(m1,m2-r2)];
   C2 = [C2 U2];
   D21 = D21+U2*So*V2';
end


% Elimination of D11 (optional)

if strcmp(option1,'D11') | strcmp(option2,'D11') | ...
      strcmp(option3,'D11'), 
   [U,S,V] = svd(D11);
   r3 = length(find(diag(S)>tol*norm([C1 D11 D12])));
   if r3 > 0
      U1 = U(:,1:r3);
      V1 = V(1:r3,:);
      S = S(1:r3,1:r3);
      n = length(E);
      E = [    E       zeros(n,r3)
            zeros(r3,n) zeros(r3,r3) ];
      A = [    A       zeros(n,r3)
         zeros(r3,n)  -eye(r3,r3) ];
      B1 = [B1; V1 ];
      B2 = [B2; zeros(r3,k2)];
      C1 = [ C1 U1*S];
      D11 = zeros(size(D11));
      C2 = [ C2 zeros(m2,r3)];
   end
end % option D11


% Elimination of D22 (optional)

if strcmp(option1,'D22') | strcmp(option2,'D22') | ...
      strcmp(option3,'D22'), 
   [U,S,V] = svd(D22);
   r3 = length(find(diag(S)>tol*norm([C2 D21 D22])));
   if r3 > 0
      U1 = U(:,1:r3);
      V1 = V(1:r3,:);
      S = S(1:r3,1:r3);
      n = length(E);
      E = [    E       zeros(n,r3)
            zeros(r3,n) zeros(r3,r3) ];
      A = [    A       zeros(n,r3)
         zeros(r3,n)  -eye(r3,r3) ];
      B1 = [B1; zeros(r3,k1) ];
      B2 = [B2; V1];
      C1 = [ C1 zeros(m1,r3) ];
      C2 = [ C2 U1*S ];
      D22 = zeros(size(D22));
   end
end % option D22


% Check

B = [B1 B2];
C = [C1; C2];
D = [D11 D12; D21 D22];

errs = [];
for omega = 1:10
   s = sqrt(-1)*omega;
   if ~(issingular(s*E-A)|issingular(s*e-a))
      G = C/(s*E-A)*B+D;
      g = c/(s*e-a)*b+d;
      errs = [errs; norm(g-G,1)/norm(g,1)];
   end
end
if verbose
   disp(sprintf('dssreg: Relative error %g',max(errs)));
elseif max(errs) > 1e-6
   disp(sprintf('dssreg warning: Relative error %g',max(errs)));
end


% Complete the operation

e = E; a = A; b = B; c = C; d = D;

%end .. dssreg
