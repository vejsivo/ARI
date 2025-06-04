function [C,u,p] = clements(P,q,tol)
%CLEMENTS  Conversion to Clements standard form
%
% The commands
%    [C,u,p] = CLEMENTS(P)
%    [C,u,p] = CLEMENTS(P,q)
%    [C,u,p] = CLEMENTS(P,q,tol)
% transform the para-Hermitian nonsingular pencil P(s) = sE+A to 
% Clements standard form
%   C(s) = u*(sE+A)*u' = se+a
% The pencil P(s) = sE+A is assumed to have no finite roots on the
% imaginary axis.
%
% The matrix u is orthogonal. The matrices e and a have the form
%    e = [ 0   0   E1
%          0   0   E3
%          E1' E3' E4 ];
%    a = [ 0   0   A1
%          0   A2  A3
%          A1' A3' A4 ];
% E1 and A1 have sizes p by p. The finite roots of the pencil 
% s*E1+A1 have strictly positive real parts. A2 is diagonal with
% the diagonal entries in ascending order.
%
% If the optional input argument q is not present then A2 has the
% largest possible size. If q is present and the largest possible
% size of A2 is greater than q by q then -- if possible -- the
% size of A2 is reduced to q by q. Setting q = Inf has the same
% effect as omitting the second input argument.
%
% The optional input parameter tol = [tol1 tol2] defines two 
% tolerances.  The tolerance tol1 is used in determining 
% whether a generalized eigenvalue of the pencil is infinite. Its 
% default value is 1e-12. The tolerance tol2 is used to test whether 
% the pencil has roots on the imaginary axis. It also has the default
% value 1e-12.
%
% If the pencil has eigenvalues on the imaginary axis then the 
% function returns p = -Inf, e = E, a = A, and u the unit matrix.
%
% Reference: D. J. Clements. "Rational spectral factorization using
% state-space methods." Systems and Control Letters, 20, pp. 335--343, 
% 1993.
%
% See also DESHINF

%    Author: H. Kwakernaak, September-October, 1997
%    Copyright(c) 1997 by Polyx, Ltd.

% Input checks

if nargin < 1,
   error('Not enough input arguments.');
end;
eval('P = pol(P);', 'error(peel(lasterr));');

if isempty(P) | size(P,1) ~= size(P,2) | P.deg ~= 1
   error('Matrix is not a square pencil.')
elseif ~isreal(P)
   error('Matrix is not real.')
elseif ~strcmp(P.v,'s') & ~strcmp(P.v,'p'),
   error('Continuous time only.')
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
   elseif ~all(size(tol)==[1 2])
      error('Invalid tolerance; must be a 1-by-2 vector.')
   end
end


% Initializations

global PGLOBAL;
eval('PGLOBAL.FORMAT;', 'painit;');
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
   tol1 = 1e-12; tol2 = 1e-12;
else
   tol1 = tol(1); tol2 = tol(2);
end


% Nonfinite deflations

nE = Inf;
while nE > length(E)
   nE = length(E);
   N = null(E);
   M = N'*A*N;
   [u,t] = schur(M);  
   Izero = find(abs(diag(t)) <= tol1*max(norm(E,1),norm(A,1)));
   if length(Izero) > 0
      v = []; nt = length(t);
      for i = Izero'
         vv = zeros(nt,1); vv(i) = 1;
         v = [v vv];
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
   end
end % Nonfinite deflations

if norm(E,1) > 1e-8*norm(A,1)
   % Finite deflation
 
   % Isolate the infinite eigenvalues
   
   [u1,e,v1] = svd(E);
   a = u1'*A*v1;
   I = find(abs(diag(e))>tol1*sum(abs(diag(e))));
   n1 = length(I); n2 = length(e)-n1;
   [u2,a2,v2] = svd(a(n1+1:n1+n2,:));
   u2 = [eye(n1) zeros(n1,n2); zeros(n2,n1) u2];
   v2 = v2(:,[n2+1:n1+n2 1:n2]);
   av2 = a*v2; a11 = av2(1:n1,1:n1);
   ev2 = e*v2; e11 = ev2(1:n1,1:n1);

   % Determine the QZ transformation
   
   [aa11,ee11,q,z] = qzord(a11,-e11,'full');
   Q = [q zeros(n1,n2); zeros(n2,n1) eye(n2)]*u2'*u1';
   Z = v1*v2*[z zeros(n1,n2); zeros(n2,n1) eye(n2)];
   
   % Separate the eigenvalues with negative real parts
   
   ad = diag(aa11); ed = diag(ee11);
   I = find(abs(ad.*conj(ed)+conj(ad).*ed)/2<tol2*abs(ed).*abs(ed));
   if length(I) > 0
      p = -Inf;  
      if verbose
         disp('clements: The pencil has roots on the imaginary axis.')
      end
      C = P; u = eye(size(P));
      return
   end
   I = find(real(diag(ee11)./diag(aa11))<-tol2);
   

   % Deflate the finite eigenvalues
   
   v = Z(:,I);
   p = length(I);
   V1 = orth([real(v) imag(v)]);
   V1 = V1(:,1:p);
   w = orth([real(Q(I,:)') imag(Q(I,:)')]);
   w = w(:,1:p);
   V2 = null([V1 w]');
   V3 = null([V1 V2]');
   V = [V1 V2 V3];
   E = V'*E*V;
   A = V'*A*V;
   nn = length(E); q = length(I);
   E = E(q+1:nn-q,q+1:nn-q);
   A = A(q+1:nn-q,q+1:nn-q);
   off = (n-length(V))/2;
   UU = eye(n);
   UU(off+1:n-off,off+1:n-off) = V;
   U = U*UU;
end


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


% Compute the final results

p = (n-length(A))/2;
C = U'*P*U;
u = U';


% Checks

Eps1 = norm(C(1:p,1:n-p),'blk',1);
Nrm = norm(P,'blk',1);
Eps = Eps1/Nrm;
if verbose
   disp(sprintf('clements: Relative residue %g',Eps));
elseif Eps > 1e-6
   disp(sprintf('clements warning: Relative residue %g',Eps));
end

%end .. clements
