function [a,b,c,d] = dss2ss(A,B,C,D,E,tol)
%DSS2SS Converts a descriptor state space system to a state space system
%
% The command
%    [a,b,c,d] = dss2ss(A,B,C,D,E [,TOL])
% converts a descriptor system to a (generalized) state
% space system. The systems are equivalent in the sense that
%
%     C (xE-A)^-1 B + D = c (xI-a)^-1 b + d(x)
%
% where x is 's', 'p' (continuous-time) or 'z', 'q' (discrete-time),
% depending on global properties setting.

%     Author: R.C.W. Strijbos, November 13, 1998.
%     Copyright 1998 by Polyx, Ltd.
%     $ Revision $  $ Date 22-Jul-2001  J.Jezek  arg check $
%                   $ Date 18-Jan-2002  J.Jezek  help      $

global PGLOBAL

eval('PGLOBAL.VARIABLE;', 'painit;');

var = PGLOBAL.VARIABLE;

switch nargin
case {5,6}
   if nargin == 5 | isempty(tol),
      tol = PGLOBAL.ZEROING;
   else
      if ~isa(tol,'double') | length(tol)>1 | ...
            ~isreal(tol) | tol<0 | tol>1,
         error('Invalid tolerance.')
      end
   end
otherwise
   error('Not enough input arguments.');
end
 
 if ~isnumeric(A) | ~isnumeric(B) | ~isnumeric(C) | ...
       ~isnumeric(D) | ~isnumeric(E) | ...
       ndims(A)>2 | ndims(B)>2 | ndims(C)>2 | ...
       ndims(D)>2 | ndims(E)>2,
    error('Invalid 1st, 2nd, 3rd, 4th or 5th argument.');
end

[rA,cA] = size(A); [rB,cB] = size(B); [rC,cC] = size(C);
[rD,cD] = size(D); [rE,cE] = size(E);

if rA~=cA | rA~=rB | rA~=rE | rA~=cC | cB~=cD | rC~=rD | rE ~= cE
   error('Matrices of inconsistent dimensions.');
end

if rA==0,
   a = A; b = B; c = C; d = D;
   return;
end;

P = shift(pol(E),1)-A;
[T,Q,Z,dims] = pencan(P,tol);
d1 = dims(1); d2 = dims(2); d3 = d1+d2;
Cz=C*Z; 
qB= Q*B;
if d1 ~= 0
   a = T(1:d1,1:d1);a=-a{0};
   b=qB(1:d1,:);
   c=Cz(:,1:d1);
else
   a=zeros(0,0);b=zeros(0,cD);c=zeros(rD,0);
end
if d2 ~= 0
   b0=qB(d1+1:d3,:);
   c0=Cz(:,d1+1:d3);
   es=eye(d2)-T(d1+1:d3,d1+1:d3);
   d=eye(d2);
   for i=1: d2-1
      d=eye(d2)+d*es;
   end
   d = c0*d*b0+D;
   d = pzer(d,tol*norm(d));
else
   d=D;
end
if isa(d,'pol')
   if isinf(d.degree) | d.degree == 0
      d = d{:};
   elseif strcmp(var,'z^-1')
      d.var = 'z';
   elseif strcmp(var,'d')
      d.var = 'q';
   end
end

%end .. dss2ss


