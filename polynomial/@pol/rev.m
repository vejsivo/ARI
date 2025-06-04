function Ar = rev(A,k,symb)
%REV   Reverse of a polynomial
%
% The command   AR = REV(A)   returns polynomial  AR(v) = v^D * A(1/v)
% where D is degree of A  and  v  is the variable of A.
% The command   AR = REV(A,K)   with scalar integer K,
% K >= DEG(A), returns polynomial  AR(v) = v^K * A(1/v) .
%
% K may also be, insteand of scalar, a matrix of the same size as A.
% For all entries, it must be  K(I,J) >= deg(A(I,J)) . In this case,
% the command acts entrywise: each K(I,J) to the corresponding A(I,J).
%
% The variable of AR is the same as that of A. If there is not any and
% someone is needed, the standard variable is used.
%
% The command   AR = REV(A,K,SYMB)  or  AR = REV(A,[],SYMB)
% acts as above but sets the variable symbol of the result to SYMB.

%       Author:  J. Jezek  17-Dec-1999
%       Copyright(c) 1999 by Polyx, Ltd.
%       $ Revision $  $ Date 06-Feb-2001 $
%                     $ Date 15-Mar-2003 $

global PGLOBAL;

na = nargin;
if na<1,
   error('Not enough input arguments.');
end;

eval('A = pol(A);','error(peel(lasterr));');
Ad = A.d;
if isempty(Ad), Ar = A; return;
end;

As = A.s; As1 = As(1); As2 = As(2);
scal = logical(1);

if na<2 | isempty(k),
   k = Ad;
else
   if ~isa(k,'double') | ~isreal(k) | ...
         any(any(floor(k)~=k | ~isfinite(k))) ,
      error('Invalid 2nd argument.');
   end;
   if length(k)==1,
      if k<Ad,
         error('Invalid 2nd argument; less than degree.');
      end;
   else
      if any(size(k)~=As),
         error('Matrices not of the same dimensions.');
      end;
      D = deg(A,'ent');
      if any(any(k<D)),
         error('Invalid 2nd argument; less than degree.');
      end;
      scal = logical(0);
   end;
end;

if na==3 & ~isempty(symb),
   if isa(symb,'pol'),
      [vs1,vs2,vd] = size(symb);
      if all([vs1,vs2,vd]==1)&(~any(symb.c(:,:)-[0,1])),
         symbh = symb.h; symb = symb.v;
      else
         error('Invalid variable symbol.');
      end;
   elseif ischar(symb),
      if strcmp(symb,'z^-1') | strcmp(symb,'zi') | ...
            (length(symb)==1 & any(symb==['z','s','p','q','d'])),
         if strcmp(symb,'zi'),
            symb = 'z^-1';
         end;
         var = pol([0 1],1,symb); symbh = var.h;
      else
         error('Invalid variable symbol.');
      end;
   else
      error('Invalid variable symbol.');
   end;   
else 
   symb = A.v; symbh = A.h;   
end; 
if isempty(symb),
   symb = PGLOBAL.VARIABLE;
   var = pol([0 1],1,symb); symbh = var.h;
end;

if isempty(Ad) | isinf(Ad),
   Ar = A; return;
end;

Ar.d = Ad;
Ar.s = As;
if scal,
   Ar.c = zeros(As1,As2,k+1);
   Ar.c(:,:,k+1-Ad:k+1) = A.c(:,:,Ad+1:-1:1);
else
   mk = max(max(k));
   for i = 1:As1,
      for j = 1:As2,
         kij = k(i,j); Dij = D(i,j);
         Ar.c(i,j,kij+1-Dij:kij+1) = A.c(i,j,Dij+1:-1:1);
      end;
   end;
end;
Ar.v = symb;
Ar.h = symbh;
Ar.u = [];
Ar.version = 3.0;

Ar = class(Ar,'pol');
Ar = pclear(Ar);

%end .. @pol/rev
