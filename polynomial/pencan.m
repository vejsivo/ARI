function [C,Q,Z,s] = pencan(P,p1,p2)
%PENCAN    Conversion to real Kronecker canonical form
%
% The function
%    [C,Q,Z,DIMS] = PENCAN(P)
%    [C,Q,Z,DIMS] = PENCAN(P[,'ord'][,TOL])
% transforms the square nonsingular real pencil P(s) = P0 + P1*s 
% to the real Kronecker canonical form 
%    C(s) = Q*P(s)*Z = diag(A+sI,I+sE)
% The real matrix A is in real Schur form, that is, it is upper 
% block triangular with the real eigenvalues on the diagonal and the 
% complex eigenvalues in 2-by-2 blocks on the diagonal.  The real 
% matrix E is upper triangular with zeros on the diagonal and, hence, 
% is nilpotent. Q and Z are real nonsingular. 
%
% If the option 'ord' is included then A is in block diagonal form 
% A = diag(A1,A2), with A1 and A2 both in real Schur form.  All
% eigenvalues of A1 have nonnegative real parts and those of A2 have 
% strictly negative real parts.
%
% The output argument DIMS contains the sizes [nA nE] of A and E. If
% the option 'ord' is included then DIMS contains the sizes [nA1 nA2 nE]
% of A1, A2 and E.
%
% The number TOL is an optional tolerance to decide which eigenvalues 
% of the pencil are infinite. Its default value is eps.
%
% If the verbose property is set to 'yes' then the macro reports the
% relative error.

% Author: H. Kwakernaak, R.C.W. Strijbos, November 13, 1998.
% Copyright 1998 by PolyX, Ltd
% Modified by J.Jezek 13-Aug-2001, arg checking

% Preparations

global PGLOBAL;
eval('PGLOBAL.VARIABLE;', 'painit;')

if nargin<1,
   error('Not enough input arguments.');
end;
eval('P = pol(P);', 'error(peel(lasterr));');
rP = size(P,1); cP = size(P,2);
if rP == cP & ~isempty(P.deg) & P.deg <= 1 
   n = cP;    
   P0 = P{0};
   if P.deg == 1
      P1 = P{1};
   else
      P1 = zeros(rP);
   end
else
   error('1st argument is not a square pencil.')
end

if rank(P) ~= rP
   error('Pencil is singular.')
end

if imag(P.coef) ~= zeros(size(P.coef)) 
   error('Pencil needs to be real.'); 
end

option = '';
epp = eps;
if nargin >= 2 
   if isstr(p1)
      option = p1;
   elseif isnumeric(p1)
      epp = p1;
   else
      error('Invalid 2nd argument.');
   end
end
if nargin == 3
   if isstr(p2)
      option = p2;
   elseif isnumeric(p2)
      epp = p2;
   else
      error('Invalid 3rd argument.');
   end
end
if ~ (strcmp(option,'ord') | strcmp(option,''))
   error('Invalid command option.'); 
end
if length(epp)~=1 | ~isreal(epp) | epp<0 | epp>1,
   error('Invalid tolerance.');
end;

ordered = strcmp(option,'ord');

% Ordered QZ transformation

if ordered
   method = 'full';
else
   method = 'partial';
end
[A,E,Q,Z] = qzord(-P0,P1,method,epp);


% Partition the pencil

I = find(abs(diag(E))>epp*max(norm(E,1),norm(A,1))); 
n1 = length(I); n2 = n-n1;

E11 = E(1:n1,1:n1); E12 = E(1:n1,n1+1:n); E22 = E(n1+1:n,n1+1:n); 
A11 = A(1:n1,1:n1); A12 = A(1:n1,n1+1:n); A22 = A(n1+1:n,n1+1:n); 
Q1 = Q(1:n1,:); Q2 = Q(n1+1:n,:); 
Z1 = Z(:,1:n1); Z2 = Z(:,n1+1:n); 


% Null the 12-block

if n1>0 & n2>0
   AA = pol([A11 E11],1); BB = pol([A22 E22],1); CC = pol(-[A12 E12],1);
   [V,W] = plyap(AA,BB,CC);
   Q1 = Q1+W*Q2;
   Z2 = Z2+Z1*V; 
end


% In the ordered case, structure the 11-block

if ordered
   d1 = diag(E); d2 = diag(A); d = d2(I)./d1(I);
   n11 = length(find(real(d)<0)); 
   n12 = length(find(real(d)>=0)); 
   e11 = E11(1:n11,1:n11); 
   e12 = E11(1:n11,n11+1:n1);
   e22 = E11(n11+1:n1,n11+1:n1);
   a11 = A11(1:n11,1:n11); 
   a12 = A11(1:n11,n11+1:n1);
   a22 = A11(n11+1:n1,n11+1:n1);
   q1 = Q1(1:n11,:); q2 = Q1(n11+1:n1,:);
   z1 = Z1(:,1:n11); z2 = Z1(:,n11+1:n1);
   if n11>0 & n12>0
      aa = pol([a11 e11],1);
      bb = pol([a22 e22],1);
      cc = pol(-[a12 e12],1);
      [v,w] = plyap(aa,bb,cc);
      E11(1:n11,n11+1:n1) = zeros(n11,n12);
      A11(1:n11,n11+1:n1) = zeros(n11,n12);
      q1 = q1+w*q2; z2 = z2+z1*v; 
   end
end


% Make E11 and -A22 equal to unit matrices

if ordered
   a11 = e11\a11; q1 = e11\q1; e11 = eye(size(e11)); 
   a22 = e22\a22; q2 = e22\q2; e22 = eye(size(e22)); 
else
   A11 = E11\A11; Q1 = E11\Q1; E11 = eye(size(E11)); 
end
E22 = -A22\E22; Q2 = -A22\Q2; A22 = -eye(size(A22));


% Make the decomposition real

if ordered
   ba = [ imag(q1); real(q1) ];
   [u,s,v] = svd(ba); uv = u(:,n11+1:2*n11)';  % left null space
   ri = uv(:,1:n11)+sqrt(-1)*uv(:,n11+1:2*n11); 
   q1 = real(ri*q1);
   a11 = real(ri*a11/ri);
   z1 = real(z1/ri);
   ba = [ imag(q2); real(q2) ];
   [u,s,v] = svd(ba); uv = u(:,n12+1:2*n12)';  % left null space
   ri = uv(:,1:n12)+sqrt(-1)*uv(:,n12+1:2*n12); 
   q2 = real(ri*q2);
   a22 = real(ri*a22/ri);
   z2 = real(z2/ri);
else
   ba = [ imag(Q1); real(Q1) ];
   [u,s,v] = svd(ba); uv = u(:,n1+1:2*n1)';  % left null space
   ri = uv(:,1:n1)+sqrt(-1)*uv(:,n1+1:2*n1); 
   Q1 = real(ri*Q1);
   A11 = real(ri*A11/ri);
   Z1 = real(Z1/ri);
end

ba = [ imag(Q2); real(Q2) ];
[u,s,v] = svd(ba); uv = u(:,n2+1:2*n2)';  % left null space
ri = uv(:,1:n2)+sqrt(-1)*uv(:,n2+1:2*n2); 
Q2 = real(ri*Q2);
E22 = real(ri*E22/ri);
Z2 = real(Z2/ri);


% Transform A11 and E33 to real Schur form

if ordered
   [u,a11] = schur(a11);
   q1 = u'*q1; z1 = z1*u;
   [u,a22] = schur(a22);
   q2 = u'*q2; z2 = z2*u;
else
   [u,A11] = schur(A11);
   Q1 = u'*Q1; Z1 = Z1*u;
end
[u,E22] = schur(E22);
Q2 = u'*Q2; Z2 = Z2*u;


% Rearrange the final result

if ordered
   Q = [q1; q2; Q2];
   Z = [z1 z2 Z2];
   AA = [     a11       zeros(n11,n12)  zeros(n11,n2)
         zeros(n12,n11)      a22        zeros(n12,n2)
         zeros(n2,n11)  zeros(n2,n12)        A22     ];
   EE = [     e11       zeros(n11,n12)  zeros(n11,n2)
         zeros(n12,n11)      e22        zeros(n12,n2)
         zeros(n2,n11)  zeros(n2,n12)        E22     ];
   s = [length(a11) length(a22) length(E22)];
else
   Q = [Q1; Q2];
   Z = [Z1  Z2];
   AA = [     A11     zeros(n1,n2)
         zeros(n2,n1)     A22     ];
   EE = [     E11     zeros(n1,n2)
         zeros(n2,n1)     E22     ];
   s = [length(A11) length(E22)];
end
C0 = -AA; C1 = EE;
C = pzer(pol([C0 C1],1),epp);

% Verbose level

verbose = strcmp(PGLOBAL.VERBOSE, 'yes');

% Checks

Eps0 = norm(Q*P*Z-C,'blk',1);
Nrm = norm(P,'blk',1);
Eps = Eps0/Nrm;
if verbose
   disp(sprintf('pencan: Relative residue %g',Eps));
elseif Eps > 1e-6
   disp(sprintf('pencan warning: Relative residue %g',Eps));
end

%end .. pencan

