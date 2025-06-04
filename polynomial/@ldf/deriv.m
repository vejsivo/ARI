function G = deriv(F,arg2,arg3,arg4)
%DERIV    Derivative of left-den fraction
%
% The command  G = DERIV(F)
% computes the derivative of left-den fraction F. The result is
% also a left-den fraction.
%
% The command  G = DERIV(F,N)
% computes the N-th derivative.
% N must be nonnegative integer sacalar.
%
% An optional input argument VAR may specify the variable for
% taking the derivative. The default is the variable of F.
% If the variable of F is 'z' or 'z^-1', VAR may be any of
% these two; otherwise, if present, VAR it must be the
% variable of F.
%
% An optional input argument TOL may specify a zeroing tolerance
% to be used instead the standard one.

%       Author:  J. Jezek  02-Feb-2000
%       Copyright(c) 2000 by Polyx, Ltd.
%       $ Revision 3.0 $  $ Date 21-Apr-2000 $
%                         $ Date 02-Nov-2000 $
%                         $ Date 30-Jul-2001  sampl per $
%                         $ Date 30-Sep-2002 $
%                         $ Date 14-Oct-2002 $
%                         $ Date 28-Feb-2003 $

global PGLOBAL;

eval('F = ldf(F);','error(peel(lasterr));');
n = 1; var = F.frac.v; h = [];
tol = PGLOBAL.ZEROING;

ni = nargin;
if ni>=2 & ~isempty(arg2),
   if isa(arg2,'char'), var = arg2;
   elseif isa(arg2,'double') & length(arg2)==1 & isreal(arg2) & arg2>=0,
      if floor(arg2)==arg2, n = arg2;
      elseif arg2<=1, tol = arg2;
      else error('Invalid 2nd argument.');
      end;
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
   if isa(arg3,'char'), var = arg3;
   elseif isa(arg3,'double') & length(arg3)==1 & isreal(arg3) & arg3>=0,
      if floor(arg3)==arg3, n = arg3;
      elseif arg3<=1, tol = arg3;
      else error('Invalid 3rd argument.');
      end;      
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
   if isa(arg4,'char'), var = arg4;
   elseif isa(arg4,'double') & length(arg4)==1 & isreal(arg4) & arg4>=0,
      if floor(arg4)==arg4, n = arg4;
      elseif arg4<=1, tol = arg4;
      else error('Invalid 4th argument.');
      end;      
   else,
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
if ~isempty(F.frac.v) & ~isempty(var) & ~strcmp(F.frac.v,var),
   if (strcmp(F.frac.v,'z') | strcmp(F.frac.v,'z^-1')) & ...
         (strcmp(var,'z') | strcmp(var,'z^-1')),
   else
      error('Invalid variable of derivative.');
   end;
end;

if ~isempty(h) & ~isempty(F.frac.h) & isfinite(h) & isfinite(F.frac.h),
   if F.frac.h~=h,
      warning('Inconsistent sampling periods.');
   end;
end;

G = F;
N = F.frac.num;
if isempty(N) | all(all(N==0)),
   return;
end;

Nder = 0; Dder = 0; X = 0; Y = 0;
if n>0,
   for i = 1:n,
      G = coprime(G);
      N = G.frac.num; D = G.frac.den;
      eval('Nder = deriv(N,var); Dder = deriv(D,var);', ...
           'error(peel(lasterr));');
      eval('[X Y] = xayb0(D,-Dder,tol);','error(peel(lasterr));');
      YNder = mtimes(Y,Nder,tol);
      XN = mtimes(X,N,tol);
      num = minus(YNder,XN,tol);
      den = mtimes(Y,D,tol);
      G = ldf(den,num);
   end;
end;

if strcmp(PGLOBAL.COPRIME,'cop'), G = coprime(G,tol);
end;
if strcmp(PGLOBAL.REDUCE,'red'), G = reduce(G,tol);
else G = smreduce(G);
end;
if strcmp(PGLOBAL.DEFRACT,'defr'), G = defract(G);
end;

%end .. @ldf/deriv

