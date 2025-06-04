function [C,u,p] = clements(P,q,tol)
%CLEMENTS1  Conversion to Clements standard form
%
% The commands
%    [C,u,p] = CLEMENTS1(P)
%    [C,u,p] = CLEMENTS1(P,q)
%    [C,u,p] = CLEMENTS1(P,q,tol)
% transform the real para-Hermitian pencil P(s) = sE+A to Clements 
% standard form
%   C(s) = u*(sE+A)*u' = se+a
% The matrix u is orthogonal. The real matrices e and a have the form
%    e = [ 0   0   E1
%          0   0   E3
%         -E1'-E3' E4 ];
%    a = [ 0   0   A1
%          0   A2  A3
%          A1' A3' A4 ];
% E1 and A1 have sizes p by p. The finite roots of the pencil 
% s*E1+A1 have nonnegative real parts. A2 is diagonal with
% the diagonal entries in ascending order.
%
% If the optional input argument q is not present then A2 has the
% largest possible size. If q is present and the largest possible
% size of A2 is greater than q by q then -- if possible -- the
% size of A2 is reduced to q by q. Setting q = Inf has the same
% effect as omitting the second input argument.
%
% The optional input parameter tol defines a relative tolerance
% with default value 1e-10. It is used to test whether 
% eigenvalues of the pencil are zero, have zero imaginary part, 
% or are infinite, and for other tests. For compatibility with an
% earlier version of the macro a tolerance parameter of the form 
% [tol1 tol2] is also accepted but only the first entry is used.
%
% The routine handles singular pencils and pencils with roots on the 
% imaginary axis.

% References: 
%
% D. J. Clements. "Rational spectral factorization using
% state-space methods." Systems and Control Letters, 20, pp. 335--343, 
% 1993.
%
% Kwakernaak, H. , "Frequency domain solution of the H-infinity
% problem for descriptor systems." In Y. Yamamoto and S. Hara, Eds.,  
% Learning, Control and Hybrid Systems, Lecture Notes in Control and 
% Information Sciences, vol. 241, Springer, London, etc., 1998.
%
% H. Kwakernaak, "A descriptor algorithm for the spectral factorization of 
% polynomial matrices," Proc. ROCOND 2000, June 21-23, 2000, Prague CZ.

%    Author: H. Kwakernaak, September-October, 1997
%    Copyright(c) 1997 by Polyx, Ltd.
%    Revised May, 1999. Revised March-April, 2000.

% Input checks

if nargin < 1,
   error('Not enough input arguments.');
end;
eval('P = pol(P);', 'error(peel(lasterr));');

if isempty(P) | size(P,1) ~= size(P,2) | P.deg~=1
   error('Matrix is not a square pencil.')
elseif ~isreal(P)
   error('Matrix is not real.')
elseif ~strcmp(P.v,'s') & ~strcmp(P.v,'p'),
   error('Continuous time only.');
elseif any(any(P'~=P))
   error('Pencil is not para-Hermitian.')
end

if nargin >= 2
   if ~isnumeric(q) | length(q)~=1 | ~isreal(q) | q<0
      error('Invalid 2nd argument.');
   end
end

if nargin == 3
   if ~isnumeric(tol) | ~isreal(tol) | ndims(tol)>2 | ...
         any(any(tol<0)) | any(any(tol>1)),
      error('Invalid 3rd argument.');
   elseif ~all(size(tol)==[1 1]) & ~all(size(tol)==[1 2]) 
      error('Invalid tolerance; must be scalar or 1-by-2 vector.')
   end
end


% Initializations

global PGLOBAL;
eval('PGLOBAL.VERBOSE;', 'painit;');
verbose = strcmp(PGLOBAL.VERBOSE, 'yes');

A = P{0}; E = P{1};
n = length(E);
U = eye(size(E));
if nargin > 1
   qq = q;
else
   qq = Inf;
end
if nargin < 3
   tol = 1e-12;
else
   tol = tol(1);
end
ep1 = tol*max(norm(E,1),norm(A,1));
ep2 = sqrt(ep1); %1e4*ep1;


% Rank deflation

[u,s,v] = svd([A;E]);
ind = find(abs(diag(s))<ep1);
V1 = v(:,ind);
V = [V1 null(V1')];
E = V'*E*V;
A = V'*A*V;  
q = size(V1,2);
nn = length(E);
E = E(q+1:nn-q,q+1:nn-q);
A = A(q+1:nn-q,q+1:nn-q);
off = (n-length(V))/2;
UU = eye(n);
UU(off+1:n-off,off+1:n-off) = V;
U = U*UU;


% Nonfinite deflations

nE = Inf;
while nE > length(E)
   nE = length(E);   
   [uu,ss,vv] = svd(E);
   II = find(diag(ss)<ep1);
   N = vv(:,II);
   M = N'*A*N;
   [u,t] = schur(M); 
   [tt,I] = sort(diag(t));
   t = t(I,I); u = u(:,I);
   tt = diag(t);
   lt = length(tt);
   if lt == 0
      break
   end
   Izero = find(abs(tt)==min(abs(tt)));
   if lt > 1 & max(abs(tt))>0
      if min(abs(tt))/max(abs(tt))>ep2, Izero = []; end
   end
   if lt == 1
      if abs(tt)>ep1, Izero = []; end
   end
   if length(Izero) > 0
      v = zeros(lt,1);
      v(Izero(1)) = 1;
   else
      break
   end
   V1 = orth(N*u*v);
   V2 = null([V1 A*V1]');
   V = [V1 V2 null([V1 V2]')];
   E = V'*E*V;
   A = V'*A*V;  
   q = size(V1,2);
   nn = length(E);
   E = E(q+1:nn-q,q+1:nn-q);
   A = A(q+1:nn-q,q+1:nn-q);
   off = (n-length(V))/2;
   UU = eye(n);
   UU(off+1:n-off,off+1:n-off) = V;
   U = U*UU;
end % Nonfinite deflations


% Deflation of the zero eigenvalues

[u,s,v] = svd(A); s = diag(s); 
while min(s) < ep1
   V1 = v(:,length(v));
   V2 = null([V1 E*V1]');
   V3 = null([V1 V2]');
   V = [V1 V2 V3];
   E = V'*E*V; A = V'*A*V;
   nn = length(E);
   E = E(2:nn-1,2:nn-1);
   A = A(2:nn-1,2:nn-1);
   off = (n-length(V))/2;
   UU = eye(n);
   UU(off+1:n-off,off+1:n-off) = V;
   U = U*UU;
   [u,s,v] = svd(A); s = diag(s);
end % Deflation of the zero eigenvalues


% Finite deflations

% Deflation of the nonzero eigenvalues on the imaginary axis

if norm(E,1) > ep2
   eigvalues = eig(A,E);
   Izero = find(abs(real(eigvalues))<ep2);
   eigzero = eigvalues(Izero);
   while length(eigzero) > 0
      lam = eigzero(1);
      [u,s,v] = svd(lam*E-A);
      I = find(diag(s) < ep2);
      e = v(:,I);
      u = real(e); v = imag(e);
      S = u'*E*v;
      [W,del] = schur(S); del = diag(del);
      [del,I] = sort(del); W = W(:,I);
      k = length(del); xi = zeros(1,k);
      if del(1) == 0 & del(k) == 0
         xi(1) = 1; xi(k) = 0;
      elseif abs(del(1)) <= abs(del(k))
         xi(1,1) = 1; xi(1,k) = sqrt(abs(del(1)/del(k)));
      else
         xi(1,1) = sqrt(abs(del(k)/del(1))); xi(1,k) = 1;
      end
      alpha = xi*W';
      u = u*alpha'; v = v*alpha';
      V1 = orth([u v]);
      if norm(V1'*E*V1,1) > ep2
         error('Cannot deflate eigenvalue on the imaginary axis.');
      end
      V2 = null([V1 E*V1]');
      V3 = null([V1 V2]');
      V = [V1 V2 V3];
      E = V'*E*V; A = V'*A*V;
      nn = length(E); q = 2;
      E = E(q+1:nn-q,q+1:nn-q);
      A = A(q+1:nn-q,q+1:nn-q);
      off = (n-length(V))/2;
      UU = eye(n);
      UU(off+1:n-off,off+1:n-off) = V;
      U = U*UU;
      eigvalues = eig(A,E);
      Izero = find(abs(real(eigvalues))<ep1);
      eigzero = eigvalues(Izero);
   end 
end % deflation of the nonzero eigenvalues on the imaginary axis


% Deflation of the finite eigenvalue off the imaginary axis

if norm(E,1) > ep2
   
   % Determine the QZ transformation
 
   [aa11,ee11,Q,Z] = qzord(A,-E,'full',ep2);

   % Separate the eigenvalues with negative real parts
   
   ad = diag(aa11); ed = diag(ee11);
   Ifinite = find(abs(ed)>ep1*abs(ad));
   eigfinite = ad(Ifinite)./ed(Ifinite); 
   Ineg = find(real(eigfinite)<-ep1);
   Izero = find(abs(real(eigfinite))<ep1);
   if length(Izero)>0
      error('Failure to extract eigenvalue on imaginary axis.')
   end
  
   
   % Deflation of the eigenvalues with negative real parts
   
   v = Z(:,Ineg);
   p = length(Ineg);
   V1 = orth([real(v) imag(v)]);
   V1 = V1(:,1:p);
   w = orth([real(Q(Ineg,:)') imag(Q(Ineg,:)')]);
   w = w(:,1:p);
   V2 = null([V1 w]');
   V3 = null([V1 V2]');
   V = [V1 V2 V3];
   E = V'*E*V;
   A = V'*A*V;
   nn = length(E); q = length(Ineg);
   E = E(q+1:nn-q,q+1:nn-q);
   A = A(q+1:nn-q,q+1:nn-q);
   off = (n-length(V))/2;
   UU = eye(n);
   UU(off+1:n-off,off+1:n-off) = V;
   U = U*UU;
end % deflation of the finite eigenvalues


if norm(E,1) < ep2

   % Make the central block diagonal
   % and make the largest entry of the corresponding
   % eigenvectors positive
 
   off = (n-length(A))/2;
   for j = off+1:n-off
      i = find(abs(U(:,j))==max(abs(U(:,j))));
      if U(i,j)<0
         U(:,j) = -U(:,j);
         A(:,j-off) = -A(:,j-off); A(j-off,:) = -A(j-off,:);
      end
   end
   [u,t] = schur(A);
   [t,I] = sort(diag(real(t))); 
   u = u(:,I); A = diag(t);
   off = (n-length(u))/2;
   UU = eye(n);
   UU(off+1:n-off,off+1:n-off) = u;
   U = U*UU;

   
   % Reduction of the central block

   while length(A) > qq 
      t = diag(A);   
      Ineg = find(t < 0); Ipos = find(t > 0);
      if min(length(Ineg),length(Ipos)) > 0
         nt = length(t);
         al = -t(1); be = t(nt);
         V = eye(nt,nt);
         V(1,1) = sqrt(be/(al+be));
         V(nt,1) = sqrt(al/(al+be));
         V(1,nt) = -V(nt,1);
         V(nt,nt) = V(1,1);
         A = V'*A*V; 
         nn = length(A);
         A = A(2:nn-1,2:nn-1);
         off = (n-length(V))/2;
         UU = eye(n);
         UU(off+1:n-off,off+1:n-off) = V;
         U = U*UU;
      else
         break
      end
   end
end % process the central block


% Compute the final results

p = (n-length(A))/2;
C = U'*P*U;
u = U';
C = pzer(C,ep1);


% Checks

Eps1 = norm(C(1:p,1:n-p),'blk',1);
Nrm = norm(P,'blk',1);
Eps = Eps1/Nrm;
if verbose
   disp(sprintf('clements: Relative residue %g',Eps));
elseif Eps > 1e-6
   disp(sprintf('clements warning: Relative residue %g',Eps));
end

%end .. clements1
