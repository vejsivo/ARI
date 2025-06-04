function [An,C] = power(A,n,tol);
%POWER(.^)   Element-wise power of polynomial
%             A.^n
%
% AN = A.^N or AN = POWER(A,N) denotes element-by-element powers. 
% The polynomial matrix A and the integer matrix N must have the
% same sizes unless one is a scalar. A scalar operates with
% every element of the other matrix.
%
% Usually, the result is also polynomial. However, if any element
% of N is negative and the corresponding element of A is non-
% constant polynomial, then the result is matrix-denominator
% fraction. The case of negative element N and corresponding
% element A equal to zero results in error. The class of the
% result can be obtained in an optional output argument.
%
% AN = POWER(A,N,TOL) works with zeroing specified by the input
% relative tolerance TOL.
%
% See also POL/MPOWER, POL/TIMES.

%       Author(s):  M. Hromcik, M. Sebek 16-2-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 09-Mar-1998 14:27:34   $
%       $Revision: 3.0 $  $Date: 11-Aug-1999 12:00:00   J. Jezek  $
%                         $Date: 21-Apr-2000 12:30:00   J. Jezek  $
%                         $Date: 07-Nov-2000 12:00:00   J. Jezek  $
%                         $Date: 19-Oct-2002            J. Jezek  $
%                         $Date: 28-Feb-2003            J. Jezek  $

% Effect on other properties:
% An.u: UserData are deleted.

global PGLOBAL;

na  = nargin;
if na<2,
   error('Not enough input arguments.');
elseif na==2 | isempty(tol),
   tol = PGLOBAL.ZEROING;
else
   if ~isa(tol,'double') | length(tol)~=1 | ...
         ~isreal(tol) | tol<0 | tol>1,
      error('Invalid tolerance.');
   end;
end;   

if ~isa(n,'double') | ndims(n)>2 | ~isreal(n) |  ...
      (~isempty(n) & any(any(floor(n)~=n)) ),
   error('Invalid power; must be integer.');
end;

degA = deg(A,'ent'); negpw = 0;
eval('negpw = any(any(n<0 & isinf(degA)));',...
     'error(peelf(lasterr));');
if negpw,
   error('Zero powered to a negative power.');
end;
if any(any(n<0 & degA>=0)),
   eval('A = mdf(A); An = power(A,n,tol);', ...
        'error(peel(lasterr));');
   C = class(An); return;
end;

me = [];

% n is scalar
if all(size(n)==1),
   An = sctimes(A,n,tol);

% n is a matrix
else
   
   if any( A.s - size(n) ),
      
     %  A is scalar
     if all(A.s==1),
       A.c = repmat(A.c, size(n));
       A.s = size(n);
     else,
       error('Matrices not of the same dimensions.');
     end;
   
   end;
   
   An = zeros(0, A.s(2) );
   Ac = A.c;
   Ad = A.d;
   Av = A.v;
   Ah = A.h;
   for i=1:A.s(1),
     an = zeros(1,0);
     for j=1:A.s(2),
       Acij = Ac(i,j,:);
       if isempty(Acij), Acij = 0; end;
       Acij = pol(Acij(:).', Ad, Av);
       Acij.h = Ah;
       an = horzcat(an,sctimes(Acij,n(i,j),tol));	
     end;
     An = vertcat(An,an);
   end;      
   
end;

if isempty(An), An = pol(An); end;
C = class(An);

% Sub-function for scalar power
function An = sctimes(A,n,tol);

nstore = n;

if n==0,
   As = A.s;
   An = pol(ones(As(1),As(2)));
% An.d = 0;
% An.s = As;
% An.c = ones(As(1),As(2),1);
% An.v = ''; 
% An.u = [];
% An.version = 3.0;
% An = class(An, 'pol');
 
elseif isinf(A.d),
 An = A;

elseif isempty(A.c),
 An = pol(zeros(A.s));
 
else 
 bin=[];
 while n~=0,
  bin= [bin, rem(n,2)];		
  n  = floor(n/2);
 end; 	%while
 Z=A;
 q=1;
 while bin(q)==0,
  Z=times(Z,Z,0);		% no zeroing yet
  q=q+1;
 end; 	%while
 F=Z;
 for k=q+1 : size(bin,2),
  Z=times(Z,Z,0);		% no zeroing yet
  if bin(k)~=0,
   F=times(F,Z,0);
  end; %if
 end; 	%for
 
 An=F;

 % zeroing
 me = min(abs(nonzeros(A.c))) ^ nstore;
 if ~isempty(me),
    An = pzer(An,tol*me);	% zeroing
 end

end;	%if else

%end .. sctimes

%end .. @pol/power
