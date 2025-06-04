function C = eq(A, B, tol);
%EQ (==)   Test if polynomials equal
%	    A == B
%
% The commmand
%    C = (A==B) 
% performs elementwise comparison between the polynomial matrices 
% A and B with a tolerance activated through the global variable
% PGLOBAL.ZEROING. A and B must have the same dimensions unless 
% one is a scalar polynomial. The scalar is compared with every
% entry of the other matrix.
%
% A difference between the polynomial variables of A and B causes a 
% warning message, but does not affect the result - for instance,
%   s == z  is 1.    
% However, if one of the variables is 'z' and the other 'z^-1',
% then the variable names play a role in the comparison;
% no warning is issued in such a case. So,
%   z == z^-1  is 0.
% 
% The commmands
%    C = EQ(A,B) 
% works alike. The commmand
%    C = EQ(A,B,TOL) 
% works with tolerance specified by the input tolerance TOL.
%  
% See also: POL/NE

%       Author(s): M. Hromcik, M. Sebek 16-2-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 21-Apr-1998 11:36:34   $
%       $Revision: 3.0 $  $Date: 11-Aug-1999 12:00:00  J.Jezek  $    
%                         $Date: 02-Aug-2000 15:00:00  J.Jezek  $
%                         $Date: 29-Jan-2002           J.Jezek  $

% Effect on other properties:
% C is a standard Matlab matrix.

global PGLOBAL;

ni = nargin;
if ni<=1,
   error('Not enough input arguments.');
elseif ni==2 | isempty(tol),
   tol = PGLOBAL.ZEROING;
else
   if ~isa(tol,'double') | length(tol)~=1 | ...
         ~isreal(tol) | tol<0 | tol>1,
      error('Invalid tolerance.');
   end;
end;

eval('A = pol(A);','error(''Invalid 1st argument.'');');
eval('B = pol(B);',...
   'eval(''B = mdf(B);'',''error(''''Invalid 2nd argument.'''');'');');
if isa(B,'mdf'),
   eval('C = eq(A,B,tol);','error(peel(lasterr));');
   return;
end;

[td,A,B] = testdp(A,B);
if td==0,
   error('Matrices not of the same dimensions.');
end;

r = zeros(A.s(1),A.s(2));
[tv,v,A,B] = testvp(A,B);
if tv == 2,
   Ad = deg(A,'ent'); Bd = deg(B,'ent');
   r = Ad>0 | Bd>0;
   A.v = 'z'; B.v = 'z';
elseif tv == 0,
   warning('Inconsistent variables.');
end;

[th,Ch,A,B] = testhp(A,B,v);
if th==0,
   warning('Inconsistent sampling periods.');
end;

   
Dif = minus(A,B,tol);
C = ~sum( cat( 3, abs(Dif.c), zeros([A.s,2]) ), 3 );
C = C & ~r;

%end .. @pol/eq
