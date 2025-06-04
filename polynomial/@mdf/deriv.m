function Q = deriv(P,arg2,arg3,arg4)
%DERIV    Derivative of matrix-den fraction
%
% The command  Q = DERIV(P)
% computes the derivative of matrix-denominator fraction P.
% The result is also a matrix-den fraction.
%
% The command  Q = DERIV(P,N)
% computes the N-th derivative.
% N must be nonnegative integer scalar.
%
% An optional input argument VAR may specify the variable for
% taking the derivative. The default is the variable of F.
% If the variable of F is 'z' or 'z^-1', VAR may be any of
% these two; otherwise, if present, VAR it must be the
% variable of F.
%
% An optional input argument TOL may specify a zeroing tolerance
% to be used instead the standard one.

%       Author:  J. Jezek  31-Jan-2000
%       Copyright(c) 2000 by Polyx, Ltd.
%       $ Revision 3.0 $   $ Date 26-Apr-2000 $
%                          $ Date 06-Nov-2000 $
%                          $ Date 30-Jul-2001  sampl per $
%                          $ Date 30_Sep-2002 $
%                          $ Date 14-Oct-2002 $
%                          $ Date 28-Feb-2003 $

global PGLOBAL;

eval('P = mdf(P);','error(peel(lasterr));');
n = 1; var = P.frac.v; h = [];
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
         var = arg2.v; h =arg2.h;
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
if ~isempty(P.frac.v) & ~isempty(var) & ~strcmp(P.frac.v,var),
   if (strcmp(P.frac.v,'z') | strcmp(P.frac.v,'z^-1')) & ...
         (strcmp(var,'z') | strcmp(var,'z^-1')),
   else
      error('Invalid variable of derivative.');
   end;
end;

if ~isempty(h) & ~isempty(P.frac.h) & isfinite(h) & isfinite(P.frac.h),
   if P.frac.h~=h,
      warning('Inconsistent sampling periods.');
   end;
end;

Q = P;
N = P.frac.num;
if isempty(N) | all(all(N==0)),
   return;
end;

Nder = 0; Dder = 0; NderD = 0;
if n>0,
   for i = 1:n,
      Q = coprime(Q,tol);
      N = Q.frac.num; D = Q.frac.den;
      eval('Nder = deriv(N,var); Dder = deriv(D,var);', ...
           'error(peel(lasterr));');
      eval('NderD = times(Nder,D,tol);','error(peel(lasterr));');
      NDder = times(N,Dder,tol);
      Q.frac.num = minus(NderD,NDder,tol);
      Q.frac.den = power(D,2,tol);
      Q.frac.c = 'cop?';
   end;
end;

if strcmp(PGLOBAL.COPRIME,'cop'), Q = coprime(Q,tol);
end;
if strcmp(PGLOBAL.REDUCE,'red'), Q = reduce(Q);
else Q = smreduce(Q);
end;
if strcmp(PGLOBAL.DEFRACT,'defr'), Q = defract(Q);
end;

%end .. @mdf/deriv
