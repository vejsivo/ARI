function Q = deriv(T,arg2,arg3,arg4)
%DERIV    Derivative of two-sided polynomial
%
% Q = DERIV(T)   computes the derivative of tsp matrix T
% Q = DERIV(T,N) computes the N-th derivative
%
% N must be nonnegative integer scalar.
%
% An optional input argument VAR  may specify the variable
% for taking the derivative. It must be 'z' or 'z^-1'.
% The default is 'z'.
%
% See also TSP/INTEGRAL.

%       Author:  J. Jezek  31-Jan-2000
%       Copyright(c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 30-Jul-2001 sampl per $
%                     $ Date 28-Feb-2003 tol $

eval('T = tsp(T);','error(peel(lasterr));');

n = 1; var = 'z'; h = [];
ni = nargin;
if ni>=2 & ~isempty(arg2),
   if isa(arg2,'double') & length(arg2)==1 & isreal(arg2) & arg2>=0,
      if floor(arg2)==arg2, n = arg2;
      elseif arg2<=1,
      else error('Invalid 2nd argument.');
      end;
   elseif isa(arg2,'char'), var = arg2;
   else
      eval('arg2 = pol(arg2);', ...
         'error(''Invalid 2nd argument.'');');      
      [vs1,vs2,vd] = size(arg2);
      if all([vs1,vs2,vd]==1) & all(arg2.c(:,:)==[0,1]),
         var = arg2.v; h = arg2.h;
      else error('Invalid 2nd argument.');
      end;
   end;
end;

if ni>=3 & ~isempty(arg3),
   if isa(arg3,'double') & length(arg3)==1 & isreal(arg3) & arg3>=0,
      if floor(arg3)==arg3, n = arg3;
      elseif arg3<=1,
      else error('Invalid 3rd argument.');
      end;
   elseif isa(arg3,'char'), var = arg3;
   else
      eval('arg3 = pol(arg3);', ...
         'error(''Invalid 3rd argument.'');');      
      [vs1,vs2,vd] = size(arg3);
      if all([vs1,vs2,vd]==1) & all(arg3.c(:,:)==[0,1]),
         var = arg3.v; h = arg3.h;
      else error('Invalid 3rd argument.');
      end;   
   end;
end;

if ni==4 & ~isempty(arg4),
   if isa(arg4,'double') & length(arg4)==1 & isreal(arg4) & arg4>=0,
      if floor(arg4)==arg4, n = arg4;
      elseif arg4<=1,
      else error('Invalid 4th argument.');
      end;
   elseif isa(arg4,'char'), var = arg4;
   else
      eval('arg4 = pol(arg4);', ...
         'error(''Invalid 4th argument.'');');
      [vs1,vs2,vd] = size(arg4);
      if all([vs1,vs2,vd]==1) & all(arg4.c(:,:)==[0,1]),
         var = arg4.v; h = arg4.h;
      else error('Invalid 4th argument.');
      end;
   end;
end;

if strcmp(var,'zi'),
   var = 'z^-1';
end;
if ~isempty(var) & ~strcmp(var,'z') & ~strcmp(var,'z^-1'),
   error('Invalid variable of derivative.');
end;

if ~isempty(h) & ~isempty(T.h) & isfinite(h) & isfinite(T.h),
   if T.h~=h,
      warning('Inconsistent sampling periods.');
   end;
end;

P = pos(T); N = neg(T);  
eval('Q = tsp(deriv(P,n,var),deriv(N,n,var));',...
   'error(peel(lasterr));');

%end .. @tsp/deriv
