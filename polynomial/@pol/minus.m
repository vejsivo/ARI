function C = minus(A,B,tol)
%MINUS  Subtract polynomials
%            C = A - B
%
% The commands
%    C = A - B    or    C = MINUS(A,B) 
% subtract the polynomial matrix B from A with zeroing using 
% the global zeroing tolerance. The commmand
%    C = MINUS(A,B,TOL) 
% works with zeroing specified by the input relative tolerance TOL. 
%
% See also POL/PLUS, POL/UMINUS.

%       Author(s):  M. Hromcik, S. Pejchova, M. Sebek 24-2-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 23-Oct-1998 10:55:34   $
%       $Revision: 3.0 $  $Date: 11-Aug-1999 12:00:00   J.Jezek  $
%                         $Date: 11-Oct-1999 12:00:00   J.Jezek  $
%                         $Date: 02-Aug-2000 11:00:00   J.Jezek  $
%                         $Date: 07-Nov-2000 12:00:00   J.Jezek  $
%                         $Date: 02-Aug-2001 12:00:00   J.Jezek  $
%                         $Date: 24-Jan-2002            J.Jezek  $

% Effect on other properties:
% C.v: Variable is kept if it is the same in both A and B, otherwise
%      it is set equal to the current value of the global variable symbol.
%      However, if one of the variables is 'z' and the other 'z^-1'
%      the result is a two-sided polynomial.
% C.u: UserData are deleted.

global PGLOBAL;

na = nargin;
if na==1,
   error('Not enough input arguments.');
elseif na==2 | isempty(tol),
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
   C = 0;
   eval('C = defract(minus(A,B,tol));','error(peel(lasterr));');
   return;
end;

[tv,Cv,A,B] = testvp(A,B);
if tv==2,
   C = 0;
   eval('C = minus(tsp(A),tsp(B),tol);','error(peel(lasterr));');
   return;
elseif tv==0,
   warning('Inconsistent variables.');
end;

[th,Ch,A,B] = testhp(A,B,Cv);
if th==0,
   warning('Inconsistent sampling periods.');
end;

[td,A,B] = testdp(A,B);
if td==0,
   error('Matrices not of the same dimensions.');
end;

me = [];
if isinf(A.d),
   C.d = B.d; C.s = B.s; C.c=-B.c;
elseif isinf(B.d),
   C.d = A.d; C.s = A.s; C.c=A.c;
else
  C.d = max(A.d, B.d);
  s = size(B.c,3) - size(A.c,3);
  C.s = A.s;
  C.c = cat(3,A.c,zeros(C.s(1),C.s(2),s)) - ...
        cat(3,B.c,zeros(C.s(1),C.s(2),-s));
  me = min(min(abs(nonzeros(A.c))),min(abs(nonzeros(B.c))));
end;

C.v = Cv;
C.h = Ch;
C.u = [];
C.version = 3.0; 
C = class(C,'pol');

% zeroing
if ~isempty(me),
   C = pzer(C,tol*me);
else
   C = pzer(C,-1);
end

%end .. @pol/minus
