function C = times(A,B,tol)
%TIMES (.*)  Element-wise multiply polynomials
%           C = A.*B
%
% C = A.*B or C = TIMES(A,B) denotes element-by-element multiplication.
% The polynomial matrices A and B must have the same sizes unless one is 
% a scalar. A scalar can be multiplied into anything. Runs with
% global zeroing tolerance.
%
% C = TIMES(A,B,TOL) works with zeroing specified by the input 
% relative tolerance TOL. 
%
% The variable symbol should be the same in both A and B. If not,
% If it is not, the variable symbol of C is set to the standard
% value. However, if one of the symbols is 'z' and the other
% 'z^-1' then the symbols play a role, the result is a two-sided
% polynomial.
%
% See also POL/MTIMES, POL/POWER.

%       Author(s): M. Hromcik, M. Sebek 16-2-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 21-Apr-1998 11:23:34   $
%       $Revision: 3.0 $  $Date: 11-Aug-1999 12:00:00   J.Jezek  $
%                         $Date: 11-Oct-1999 12:00:00   J.Jezek  $
%                         $Date: 31-Jan-2000  M.Hromcik, J.Jezek $
%                         $Date: 12-Jul-2000  J.Jezek  $
%                         $Date: 02-Aug-2000  J.Jezek  $
%                         $Date: 07-Nov-2000  J.Jezek  $
%                         $Date: 24-Jan-2002  J.Jezek  $
%                         $Date: 28-Feb-2003  J.Jezek  $

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
   eval('C = defract(times(A,B,tol));','error(peel(lasterr));');
   return;
end;

me = 0;
Ad = A.d; Bd = B.d;

[tv,Cv,A,B] = testvp(A,B);
if tv==2,
   C = 0;
   eval('C = times(tsp(A),tsp(B),tol);','error(peel(lasterr));');
   return;
elseif tv==0,
   warning('Inconsistent variables.');
end;   

[th,Ch,A,B] = testhp(A,B,Cv);
if th==0,
   warning('Inconsistent sampling periods.');
end;

isscalarA = 0;
isscalarB = 0;
if any(A.s - B.s),					% fast for scalar multiple
   if all(A.s==1),
     isscalarA = 1;
   elseif all(B.s==1),
     isscalarB = 1;
   else,
     error('Matrices not of the same dimensions.');
   end;
end;

if isempty(A.d) |  isempty(B.d),       % result is empty
   C.d = [];
   if ~isscalarA,
      C.s = A.s; C.c = zeros(A.s(1),A.s(2),0);
   else
      C.s = B.s; C.c = zeros(B.s(1),B.s(2),0);
   end;
   
elseif isinf(A.d),					   	% A is zero
   C.d = -Inf;
   if ~isscalarA,
      C.s = A.s; C.c = A.c;
   else
      C.s = B.s; C.c = zeros(B.s(1),B.s(2),0);
   end;
   
elseif isinf(B.d),		      			% B is zero
   C.d = -Inf;
   if ~isscalarB,
      C.s = B.s; C.c = B.c;
   else
      C.s = A.s; C.c = zeros(A.s(1),A.s(2),0);
   end;
   
elseif (A.d==0) & (B.d==0),
   C.d = 0;
   C.s = size(A.c .* B.c);
   C.c = A.c .* B.c;
   
elseif isscalarA,
      C.d = A.d+B.d;
      C.s = B.s;
      B.c = cat(3, flipdim(B.c,3), zeros([B.s, A.d]) );
      Cc = filter( fliplr(A.c(:,:)), 1,  B.c, [], 3);
      C.c = flipdim(Cc,3);

elseif isscalarB,
      C.d = A.d + B.d;
      C.s = A.s;
      A.c = cat(3, flipdim(A.c,3), zeros([A.s, B.d]) );
      Cc = filter( fliplr(B.c(:,:)), 1,  A.c, [], 3);
      C.c = flipdim(Cc,3);

else
   C.d=A.d+B.d;
   C.s=A.s;
   if A.d * B.d < A.s(1) * A.s(2),			% large size, low degree
    C.c=zeros([C.s,C.d+1]);
    Apom=zeros([A.s, B.d+1]);

    for i=1 : A.d+1,
     Apom=repmat(A.c(:,:,i),[1 1 B.d+1]);			
     C.c(:,:,i:i+B.d)=C.c(:,:,i:i+B.d)+Apom.*B.c(:,:,:);
    end;
 
   else							%small size, high degree
    siz=A.s(1)*A.s(2);
    pom=zeros([siz, C.d+1]);
    Ar = reshape(A.c,[siz , A.d+1]);
    Br = reshape(B.c,[siz , B.d + 1]);
    for i=1:siz,
     pom(i,:) = conv( Ar(i,:),Br(i,:) );
    end; 
    C.c = reshape(pom, [C.s, C.d+1]);
   
   end; 

end;   

C.v = Cv;
C.h = Ch;
C.u = [];
C.version = 3.0;  
C = class(C,'pol');

% zeroing
if ~isempty(A.c) & ~isempty(B.c),
   me = min(min(abs(nonzeros(A.c))),min(abs(nonzeros(B.c))));
end   
C = pzer(C,tol*me);

%end .. @pol/times
