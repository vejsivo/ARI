function r = ismonomod(T,tol)
%ISMONOMOD  Test if two-sided polynomial is monomodular
%
% The command
%    R = ISMONOMOD(T[,TOL])
% returns 1 if the square tsp matrix T is monomodular and 0 if it is not.
%
% The monomodularity of tsp means: the determinant is invertible, its
% inverse being also tsp, i.e. the determinant is monomial  K*z^k
% with K nonzero, k any integer.
%
% An optional tolerance TOL may be included. Its default value is the 
% global zeroing tolerance.
%
% See also TSP/ISUNIMOD, POL/ISMONOMOD.

%       Author:  J. Jezek  30-Sep-2002
%       Copyright (c) 2002 by Polyx, Ltd.
%       $ Revision $  $ Date 28-Feb-2003  warning  $

global PGLOBAL;

if nargin<2 | isempty(tol),
   tol = PGLOBAL.ZEROING;
else
   if ~isa(tol,'double') | length(tol)~=1 | ...
      ~isreal(tol) | tol<0 | tol>1,
      error('Invalid tolerance.');
   end;
end;

if ~isa(T,'tsp'),
   error('Some argument but not 1st is invalidly tsp.');
end
[m,n] = size(T);
if m~=n,
   error('Matrix is not square.');
end;

% Test monomodularity
P = T.p; detP = det(P,0);
detPc = detP.c; detPd = detP.d;
if ne(detP,0,tol) & all(isfinite(detPc(:)))
   % det P is nonzero with tolerance TOL
   [pm,k] = max(abs(detPc(:)));
   aPc1 = abs(detPc(1:k-1)); aPc2 = abs(detPc(k+1:detPd+1));
   tolm = tol*pm;
   r = all(aPc1<tolm) & all(aPc2<tolm);
   return;
end;

% First try the rank using numerically reliable method:
rankP = rank(P,tol);
if rankP<min(m,n),  % Singular matrix
   r = 0;
   warning('Matrix is singular to working precision.');
   return;
end;
   
% P is badly scaled - is non singular and has zero determinant
% Try to scale the input P to make its determinant 
% "reasonable" - with unity leading or closing coefficient:
  
Plc = lcoef(P,'col');
saved_w = warning; warning off;
Plci = inv(Plc); warning(saved_w);
if all(all(isfinite(Plci))),     % col reduced
   P = P*Plci; detP = det(P,0);
   detPc = detP.c; detPd = detP.d;
   [pm,k] = max(abs(detPc(:)));
   aPc1 = abs(detPc(1:k-1)); aPc2 = abs(detPc(k+1:detPd+1));
   tolm = tol*pm;
   r = all(aPc1<tolm) & all(aPc2<tolm);
   return;
end;
   
Pcc = tcoef(P,'col');
savd_w = warning; warning off;
Pcci = inv(Pcc); warning(saved_w);
if all(all(isfinite(Pcci))),     % back col reduced
   P = P*Pcci; detP = det(P,0);
   detPc = detP.c; detPd = detP.d;
   [pm,k] = max(abs(detPc(:)));
   aPc1 = abs(detPc(1:k-1)); aPc2 = abs(detPc(k+1:detPd+1));
   tolm = tol*pm;
   r = all(aPc1<tolm) & all(aPc2<tolm);
   return;
end;
   
Plr = lcoef(P,'row');
saved_w = warning; warning off;
Plri = inv(Plr); warning(saved_w);;
if all(all(isfinite(Plri))),     % row reduced
   P = P*Plri; detP = det(P,0);
   detPc = detP.c; detPd = detP.d;
   [pm,k] = max(abs(detPc(:)));
   aPc1 = abs(detPc(1:k-1)); aPc2 = abs(detPc(k+1:detPd+1));
   tolm = tol*pm;
   r = all(aPc1<tolm) & all(aPc2<tolm);
   return;
end;
   
Pcr = tcoef(P,'row');
saved_w = warning; warning off;
Pcri = inv(Pcr); warning(saved_w);
if all(all(isfinite(Pcri))),     % back row reduced
   P = P*Pcri; detP = det(P,0);
   detPc = detP.c; detPd = detP.d;
   [pm,k] = max(abs(detPc(:)));
   aPc1 = abs(detPc(1:k-1)); aPc2 = abs(detPc(k+1:detPd+1));
   tolm = tol*pm;
   r = all(aPc1<tolm) & all(aPc2<tolm);
   return;
end;
  
r = 0;
warning('Matrix is badly scaled.');

%end .. @tsp/ismonomod
