function C = mtimes(A,B,tol)
%MTIMES (*)  Matrix multiply polynomials
%           C = A*B
%
% C = A*B or C = MTIMES(A,B) is the matrix product of the polynomial
% matrices A and B. Any scalar polynomial may multiply anything. 
% Otherwise, the number of columns of A must equal the number of
% rows of B. Runs with zeroing using global zeroing tolerance. 
%
% C = MTIMES(A,B,TOL) works with zeroing specified by the input 
% relative tolerance TOL. 
%
% See also POL/TIMES.

%       Author(s):  D. Henrion, S. Pejchova, M. Sebek 02-3-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 24-Sep-1998 12:12:34   $
%       $Revision: 3.0 $  $Date: 11-Aug-1999 12:00:00   J.Jezek  $
%                         $Date: 11-Oct-1999 12:00:00   J.Jezek  $
%                         $Date: 02-Aug-2000 15:00:00   J.Jezek  $
%                         $Date: 07-Nov-2000 12:00:00   J.Jezek  $
%                         $Date: 24-Jan-2002            J.Jezek  $
%       Modified by M. Sebek, 03-Dec-2011:  KRON replaced by PKRON

% Effect on other polynomial properties:
% C.v: Variable is kept if it is the same in both A and B, otherwise
%      it is set equal to the current value of the global variable symbol.
%      However, if one of the variables is 'z' and the other 'z^-1',
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
   'eval(''B = sdf(B);'',''error(''''Invalid 2nd argument.'''');'');');
if isa(B,'sdf'),
   C = 0;
   eval('C = defract(mtimes(A,B,tol));','error(peel(lasterr));');
   return;
end;

[tv,Cv,A,B] = testvp(A,B);
if tv==2,
   C = 0;
   eval('C = mtimes(tsp(A),tsp(B),tol);','error(peel(lasterr));');
   return;
elseif tv==0,
   warning('Inconsistent variables.');
end;

[th,Ch,A,B] = testhp(A,B,Cv);
if th==0,
   warning('Inconsistent sampling periods.');
end;

me = 0;
tr = 1; RA = []; RB = [];

if A.d > 1.2*(B.d),  tr = 0;
end;

C.d = A.d + B.d;

if (isempty(C.d))|(isinf(C.d)), tr = 2;
end;

if (A.s(2) ~= B.s(1)),
   if all(A.s==1),
      if tr==1, RA = pkron(A.c, eye(B.s(1))); % MS: KRON replaced by PKRON
         % build block row matrix [A0 A1 .. An]
      elseif tr==0,
         A.c = reshape(pkron(A.c, eye(B.s(1))),[B.s(1), B.s(1), A.d+1]);
         % MS: KRON replaced by PKRON
      end;
      A.s = [B.s(1) B.s(1)];
   elseif all(B.s==1),
      if tr==1,
         B.c = reshape(pkron(B.c, eye(A.s(2))),[A.s(2), A.s(2),max(0,B.d)+1]);
         % MS: KRON replaced by PKRON
      elseif tr==0,
         RB = pkron(B.c(:), eye(A.s(2))); % MS: KRON replaced by PKRON
      end;
      B.s = [A.s(2) A.s(2)];
   else,
      error('Matrices of inconsistent dimensions.');
   end;
end;

C.s = [A.s(1) B.s(2)];

if tr==1,
   if isempty(RA), RA = A.c(:,:);
   end;
   RB = sylv(B,A.d);
   RC = RA * RB;
   C.c = reshape(RC,[C.s(1),C.s(2),C.d+1]);
elseif tr==0,
   RA = sylv(A,B.d,'col');
   if isempty(RB),
      RB = reshape(permute(B.c,[1 3 2]),[(max(0,B.d)+1)*B.s(1), B.s(2)]);
   end;
   RC = RA * RB;
   C.c = permute(reshape(RC,[C.s(1),C.d+1,C.s(2)]),[1,3,2]);
else,
   if all(C.s),  C.d = -Inf;
   end;
   C.c = zeros(C.s(1), C.s(2), max(C.d+1,0));
end;

if tr~=2,
   me = min(abs(nonzeros(A.c)))*(min(abs(nonzeros(B.c))));
end;

C.v = Cv;
C.h = Ch;
C.u = [];
C.version = 3.0; 
C = class(C,'pol');

% zeroing
C = pzer(C,tol*me);

%end .. @pol/mtimes
