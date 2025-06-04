function C = kron(A,B,tol)
%KRON   Kronecker tensor product of polynomials
%
% KRON(X,Y) is the Kronecker tensor product of X and Y.
% The result is a large matrix formed by taking all possible
% products between the elements of X and those of Y.  
%
% KRON(X,Y,TOL) works with zeroing specified by the input 
% relative tolerance TOL.
%
% The variable symbol of X and Y should be the same. If not,
% it is set to the current value of the global variable
% symbol. However, if one symbol is 'z' and the other 'z^-1'
% then the symbols play a role and the result is two-sided
% polynomial.

% See also POL/TIMES, POL/MTIMES.

%       Author(s): M. Hromcik, M. Sebek 16-2-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 17-Apr-1998 14:30:34   $
%       $Revision: 3.0 $  $Date: 11-Aug-1999 12:00:00  J.Jezek  $
%                         $Date: 11-Oct-1999 12:00:00  J.Jezek  $
%                         $Date: 12-Jul-2000 14:30:00  J.Jezek  $
%                         $Date: 02-Aug-2000 17:00:00  J.Jezek  $
%                         $Date: 02-Jul-2001 J.Jezek, M.Hromcik $
%                                        $ Bug in M6.0's filter $
%                         $Date: 25-Jan-2002 J.Jezek            $

% Effect on other properties:
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
   eval('C = defract(kron(A,B,tol));','error(peel(lasterr));');
   return;
end;

[tv,Cv,A,B] = testvp(A,B);
if tv==2,
   eval('C = kron(tsp(A),tsp(B),tol);','error(peel(lasterr));');
   return;
elseif tv==0,
   warning('Inconsistent variables.');
end;

[th,Ch,A,B] = testhp(A,B,Cv);
if th==0,
   warning('Inconsistent sampling periods.');
end;

me = 0;
Ad = A.d;
Bd = B.d;
As1 = A.s(1);
As2 = A.s(2);
Bs1 = B.s(1);
Bs2 = B.s(2);

if isempty(Ad) | isempty(Bd) | Ad<0 | Bd<0,	% A or B is zero or empty
  if isempty(Ad) | isempty(Bd),
    C.d = [];
  else 
    C.d = -Inf;
  end;
  C.s = B.s .* A.s; C.c = zeros(C.s(1),C.s(2),0);

else
  C.d = A.d+B.d;
  C.s = B.s .* A.s;
  Cc = zeros( [B.s .* A.s, C.d+1] );
  Ac = permute(A.c, [3,1,2]);
  Ac = Ac(:,:); 
  Ac = flipud( Ac );
  Aci = zeros(1,Ad+1);
  B.c = cat(3, flipdim(B.c,3), zeros([B.s,A.d]) );
  for i = 1:prod(A.s),
      Aci = Ac(:,i,:); 
      vert = floor( (i-1)/As1);
      horz = rem( (i-1),As1 ); 
      %eval('Cc( horz*Bs1+1:(horz+1)*Bs1, vert*Bs2+1:(vert+1)*Bs2,:) = filter( Aci, 1,  B.c, [], 3);', ...
      %     'Cc( horz*Bs1+1:(horz+1)*Bs1, vert*Bs2+1:(vert+1)*Bs2,:) = Aci*B.c;');
      Bcpom = permute(B.c,[3,1,2]);
      pom = filter( Aci, 1,  Bcpom);
      Cc( horz*Bs1+1:(horz+1)*Bs1, vert*Bs2+1:(vert+1)*Bs2,:) = permute(pom,[2 3 1]);      
  end;
  C.c = flipdim( reshape(Cc, [C.s, C.d+1]), 3 );
end;  

C.v = Cv;
C.h = Ch;
C.u = [];
C.version = 3.0;  
C = class(C,'pol');

% zeroing
if ~isempty(A.c) & ~isempty(B.c),
   me = min(min(abs(nonzeros(A.c))),min(abs(nonzeros(B.c))));
end;
C = pzer(C,tol*me);

%end .. @pol/kron
