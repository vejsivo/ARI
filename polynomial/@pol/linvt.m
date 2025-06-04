function Q = linvt(P, a, b);
% LINVT  Linear transformation of the variable of a polynomial matrix
%
% If P is a polynomial matrix with variable VAR and A and B are real 
% numbers, then
%    Q = LINVT(P,A,B)
% computes the polynomial matrix Q such that 
%	 Q(VAR) = P(A*VAR + B)
%
% See also POL/SHIFT, POL/REVERSE, POL/SUBSREF.

%       Author(s): M. Hromcik, M. Sebek 3-9-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 03-Sep-1998 10:28:34   $
%       $Revision: 3.0 $  $Date: 30-May-2000  J.Jezek   $ 

% Algorithm: cf. D. Bini, V. Pan: Polynomial and Matrix Computations,
%	     Volume 1, Fundamental Algorithms, Birkhauser, Boston 1994,
%	     pp. 15, 16.

if nargin == 0,
   error('Not enough input arguments.');
end;

eval('P = pol(P);', 'error(peel(lasterr));');

if nargin == 1,
  a = 1; b = 0;
elseif nargin == 2,
  b = 0;
end;    

if ~isa(a, 'double') | ~isa(b, 'double') | ...
    any( size(a) > 1 ) | any( size (b) >1 ),
  
  error('Invalid 2nd or 3rd argument; must be scalars.');

end;

n = P.d;
Pc = P.c;
Ps = P.s;

if isempty(n) | n <= 0,		% zero & empty & constant
  Q = P;
  return;
end;    

u = zeros([Ps, n+1]);
v = zeros(1, n+1);
fact = zeros(1, n+1);

for h = 0:n,
  fact_h = prod(1:h);
  u(:,:,n-h+1) = fact_h * Pc(:,:,h+1);
  v(h+1) = b^h / fact_h;
  fact(h+1) = fact_h;
end;

%w .. conv(u, v):
u(:,:,length(v)+size(u,3)-1) = zeros(Ps);
w = filter(v,1,u,[],3);

for g = 0:n,
  q(:,:,g+1) = a^g * w(:,:,n-g+1) / fact(g+1);
end;  

Q.d = n; Q.s = Ps; Q.c = q; Q.v = P.v; Q.h = P.h; 
Q.u = []; Q.version = 3.0;
Q = class(Q, 'pol');

%end .. @pol/linvt
