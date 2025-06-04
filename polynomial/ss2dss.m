function [a,b,c,d,e] = ss2dss(A,B,C,D,tol)
%SS2DSS  Convert a state space system to a descriptor state space system
%
% The command
%    [a,b,c,d,e] = ss2dss(A,B,C,D [,TOL])
% with D a constant matrix or a polynomial matrix, converts the (generalized) 
% state space representation to a descriptor state space representation.
% The systems are equivalent in the sense that
%     C (xI-A)^-1 B + D(x) = c (xe-a)^-1 b + d
% where x is 's' or 'p' (continuous-time) or 'z' or 'q' (discrete-time).
% Here a,b,c,d and e are constant matrices.

%     Author: R.C.W. Strijbos, November 13, 1998.
%     Copyright 1998 by Polyx, Ltd.
%     $ Revision $  $ Date 22-Jul-2001  J.Jezek  arg check $
%                   $ Date 18-Jan-2002  J.Jezek  help      $


switch nargin
case {4,5}
   if nargin == 4
      tol = eps*10^2;
   else
      if ~(isa(tol,'double') & length(tol)==1 & ...
            isreal(tol) & tol>=0 & tol<=1),
         error('Invalid tolerance.')
      end
   end
otherwise
   error('Not enough input arguments.');
end

if ~isnumeric(A) | ~isnumeric(B) | ~isnumeric(C) | ...
      ndims(A)>2 | ndims(B)>2 | ndims(C)>2,
   error('Invalid 1st, 2nd or 3rd argument.');
end;

eval('D = pol(D);', ...
   'error(''Invalid 4th argument; not convertible to polynomial.'');');
var = D.v;
if strcmp(var,'z^-1') | strcmp(var,'d'),
   error('Invalid variable symbol in 4th argument.');
end;

[rA,cA] = size(A); [rB,cB] = size(B);
[rC,cC] = size(C); [rD,cD] = size(D);
if rA~=cA | rA~=rB | rA~=cC | cB~=cD | rC~=rD
   error('Matrices of inconsistent dimensions.');
end

d = D{0};
T = [eye(rD) shift(D-d,1)];
T = rowred(T,tol);
DegT = deg(T,'row');
for i = 1:rD
   k = T(i,:);
   k = flipud(reshape(k{:},rD+cD,DegT(i)+1)')';
   T(i,:) = pol(k(:)',DegT(i),T.var);
end
nh = T(:,rD+1:rD+cD);
dh = rowred(T(:,1:rD));
[a1,b1,c1,d1] = lmf2ss(nh,dh);
la = length(a1);
a = A;b = B;c = C;
a(rA+1:rA+la,rA+1:rA+la) = eye(la);
b=[B;b1];
c=[C -c1];
e=eye(rA);
e(rA+1:rA+la,rA+1:rA+la) = a1;

%end .. ss2dss

