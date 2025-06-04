function [M,C] = inv(T,arg2,arg3)
%INV    Inverse of two-sided polynomial
%
% The command  M = INV(T)  returns the inverse of square two-sided
% polynomial matrix T. If possible, the result is also a two-sided
% polynomial, if not, it is a scalar-denominator fraction in 
% variable 'z' or 'z^-1', decided standardly by PGLOBAL.DISCRVAR
% or by optional input argument VAR. The class of the result can
% be obtained in an optional output argument.
%
% An optional input argument TOL may specify a zeroing tolerance
% to be used instead of PGLOBAL.ZEROING.

%        Author:  J. Jezek  24-Mar-2000
%        Copyright(c) 2000 by Polyx, Ltd.
%        $ Revision 3.0 $  $ Date 21-Apr-2000 $
%                          $ Date 23-May-2000 $
%                          $ Date 31-Oct-2000 $

global PGLOBAL;

tol = []; var = '';
if nargin>=2,
   if isa(arg2,'double'), tol = arg2;
   elseif isa(arg2,'char'), var = arg2;
   elseif isa(arg2,'pol'),
      [vs1,vs2,vd] = size(arg2);
      if all([vs1,vs2,vd]==1) & all(arg2.c(:,:)==[0,1]),
         var = arg2.v;
      else error('Invalid 2nd argument.');
      end;
   else error('Invalid 2nd argument.');
   end;
   if nargin==3,
      if isa(arg3,'double'), tol = arg3;
      elseif isa(arg3,'char'), var = arg3;
      elseif isa(arg3,'pol'),
         [vs1,vs2,vd] = size(arg3);
         if all([vs1,vs2,vd]==1) & all(arg3.c(:,:)==[0 1]),
            var = arg3.v;
         else error('Invalid 3rd argument.');
         end;
      else error('Invalid 3rd argument.');
      end;
   end;
end;

if isempty(tol),
   tol = PGLOBAL.ZEROING;
elseif length(tol)~=1 | ~isreal(tol) | tol<0 | tol>1,
   error('Invalid tolerance.');
end;

if isempty(var),
   var = PGLOBAL.DISCRVAR;
elseif ~strcmp(var,'z') & ~strcmp(var,'zi') & ~strcmp(var,'z^-1'),
   error('Invalid variable symbol.');
end;

if isempty(T),
   if size(T,1)~=size(T,2),
      error('Matrix is not square.');
   end;
   M = T; C = class(M); return;
end;

[dT,lc]= deg(T); M = 0;
if tdeg(T)==dT,
   eval('M = tsp(inv(lc)*z^-dT);','error(peelf(lasterr));');
   M.h = T.h; M = tclear(M);
   C = class(M); return;
end;
   
eval('M = inv(T.p,tol);','error(peel(lasterr));');

[dD,lcD] = deg(M.d);
if tdeg(M.d)==dD,
   M = (shift(tsp(M.n),-(dD+T.o)))/lcD;
   M.h = T.h; M = tclear(M);
   C = class(M); return;
end;

if T.o>=0, M.d = shift(M.d,T.o,'z');
else       M.n = shift(M.n,abs(T.o),'z');
end;
if ~strcmp(var,'z'), M = reverse(M);
end;
C = class(M);

%end .. @tsp/inv

