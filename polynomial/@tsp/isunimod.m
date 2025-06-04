function r = isunimod(T,arg2)
%ISUNIMOD Test if two-sided polynomial matrix is unimodular
%
% The command
%    R = ISUNIMOD(T[,TOL])
% returns 1 if the two-sided polynomial matrix T is unimodular
% and 0 if it is not.
%
% An optional tolerance TOL may be included. Its default value is the 
% global zeroing tolerance.
%
% See also POL/ISUNIMOD, TSP/ISMONOMOD.

%	Author(s): M. Hromcik, H. Kwakernaak, M. Sebek 30-10-98,  J. Jezek
%	Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 19-Nov-1998 10:28:34   $
%       $Revision: 3.0 $  $Date: 28-Feb-2003  J.Jezek  warning  $

global PGLOBAL;

eval('tol = PGLOBAL.ZEROING;', 'painit; tol = PGLOBAL.ZEROING;'); 

switch nargin,
case 0, error('Not enough input arguments.');
case 1,
case 2,
   if ~isempty(arg2),
      if isnumeric(arg2) & length(arg2)==1 & isreal(arg2) & ...
            arg2>=0 & arg2 <=1,
         tol = arg2;
      else
         error('Invalid tolerance.');
      end;
   end;    
end;		%switch

% Initializations
if ~isa(T,'tsp'),
   error('Some argument but not 1st is invalidly tsp.');
end;
[m,n] = size(T);
if n ~= m
   error('Matrix is not square.')
end;
   
% Test unimodularity:
P = T.p; detP = det(P,0); mTo = m*T.o;
detPc = detP.c; detPd = detP.d;
if ne(detP,0,tol) & all(isfinite(detPc(:))) 
   % det P is nonzero with tolerance TOL
   [pm,k] = max(abs(detPc(:)));
   if mTo + k-1 ~=0,
      r = logical(0); return;
   end;
   aPc1 = abs(detPc(1:k-1)); aPc2 = abs(detPc(k+1:detPd+1));
   tolm = tol*pm;
   r = all(aPc1<tolm) & all(aPc2<tolm);
   return;
end;

% First try the rank using numerically reliable method:
switch nargin,
case 1,
  	rankP = rank(P);
case 2,
	rankP = rank(P,tol);
end;		%switch

if rankP < min(m,n),	% Singular matrix
   r = 0;
   warning('Matrix is singular to working precision.');
   return;
end;  
  
% P is badly scaled - is non singular and has zero determinant
% Try to scale the input P to make its determinant 
% "reasonable" - with unity leading or closing coefficient   
    
Plc = lcoef(P,'col');
saved_w = warning; warning off;
Plci = inv(Plc); warning(saved_w);
if all(all(isfinite(Plci))),     % col reduced
   P = P*Plci; detP = det(P,0);
   detPc = detP.c; detPd = detP.d;
   [pm,k] = max(abs(detPc(:)));
   if mTo + k-1 ~=0,
      r = logical(0); return;
   end;
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
   if mTo + k-1 ~=0,
      r = logical(0); return;
   end;   
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
   if mTo + k-1 ~=0,
      r = logical(0); return;
   end;   
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
   if mTo + k-1 ~=0,
      r = logical(0); return;
   end;   
   aPc1 = abs(detPc(1:k-1)); aPc2 = abs(detPc(k+1:detPd+1));
   tolm = tol*pm;
   r = all(aPc1<tolm) & all(aPc2<tolm);
   return;
end;

r = 0;
warning('Matrix is badly scaled.');

%end .. @tsp/isunimod
 
