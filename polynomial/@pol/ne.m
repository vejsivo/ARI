function C = ne(A, B, tol);
%NE (~=)   Test if polynomials unequal
%	     A ~= B	
%
% C = A~=B does elementwise comparison of the polynomial matrices 
% A and B with global zeroing tolerance. A and B must have the same 
% dimensions unless one is a scalar polynomial. A scalar can be 
% compared with anything.
%
% A difference between the polynomial variables of A and B causes
% a warning  message, but does not affect the result - for instance,
%   z ~= s  is 0.  
% However, if one of the variables is 'z' and the other 'z^-1',
% then the variable names play a role in the comparison;
% no warning is issued in such a case. So,
%   z ~= z^-1  is 1.
% 
% C = NE(A,B) works alike.
% C = NE(A,B,TOL) works with tolerance specified by the input tolerance TOL.
% 
%    See also: POL/EQ.

%       Author(s): M. Hromcik, M. Sebek 16-2-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 21-Apr-1998 11:41:34   $
%       $Revision: 3.0 $  $Date: 11-Aug-1999 12:00:00   J.Jezek  $
%                         $Date: 02-Aug-2000 15:00:00   J.Jezek  $
%                         $Date: 29-Jan-2002            J.Jezek  $

% Effect on other properties:
% C is a standard Matlab matrix.

global PGLOBAL;

ni = nargin;
if ni<=1,
   error('Not enough input arguments.');
elseif ni==2  | isempty(tol),
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
   eval('C = ne(A,B,tol);','error(peel(lasterr));');
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

[th,Ch,A,B] = testhp(A,B,'');
if th==0,
   warning('Inconsistent sampling periods.');
end;

Dif = minus(A,B,tol);
C = 1 & sum( cat( 3, abs(Dif.c), zeros([A.s,2]) ), 3);
C = C | r;

%end .. @pol/ne


 

