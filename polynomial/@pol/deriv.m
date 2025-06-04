function Q = deriv(P,arg2,arg3,arg4)
%DERIV    Derivative of polynomial
%
% Q = DERIV(P)   computes the derivative of polynomial matrix P
% Q = DERIV(P,N) computes the N-th derivative
%
% N must be nonnegative integer scalar.
%
% An optional input argument VAR may specify the variable
% for taking the derivative. The default is the variable
% of P. If P is polynomial in 'z' or in 'z^-1', VAR may be
% any of these two variables; otherwise, if present, it
% must be the variable of P.
%
% See also TSP/DERIV, POL/INTEGRAL.

%       Author:  J. Jezek  31-Jan-2000
%       Copyright(c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 09-Oct-2000 $
%                     $ Date 30-Jul-2001 sampl per $
%                     $ Date 28-Feb-2003 tol $

eval('P = pol(P);','error(peel(lasterr));');
n = 1; var = P.v; h = [];

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
if ~isempty(P.v) & ~isempty(var) & ~strcmp(P.v,var),
   if (strcmp(P.v,'z') | strcmp(P.v,'z^-1')) & ...
         (strcmp(var,'z') | strcmp(var,'z^-1')),
   else
      error('Invalid variable of derivative.');
   end;
end;

if ~isempty(h) & ~isempty(P.h) & isfinite(h) & isfinite(P.h),
   if P.h~=h,
      warning('Inconsistent sampling periods.');
   end;
end;

if isempty(P),
   Q = P; return;
end;

Pd = P.d; Ps = P.s; Ps1 = Ps(1); Ps2 = Ps(2);
Q = pol(zeros(Ps));

if strcmp(P.v,var) | isempty(P.v) | isempty(var),
   if n<=Pd,
      numbs = 1:Pd;
      Pd1 = Pd+1; coefs = ones(1,Pd1);
      if n>0,
         for i = 1:n,
            coefs = coefs(2:end).*numbs(1:Pd1-i);
         end;
      end;
      coefs = kron(coefs,ones(Ps));
      coefs = reshape(coefs,Ps(1),Ps(2),Pd1-n);
      Q.c = P.c(:,:,1+n:Pd1).*coefs;
      Q.d = Pd-n;
   end;
   
else
   Pdm1 = Pd-1;
   numbs = -(1:Pdm1+n);
   coefs = ones(1,Pd);
   if n>0,
      for i = 1:n,
         coefs = coefs.*numbs(i:i+Pdm1);
      end;
   end;
   coefs = kron(coefs,ones(Ps));
   coefs = reshape(coefs,Ps1,Ps2,Pd);
   Q.c = cat(3,zeros(Ps1,Ps2,1+n),P.c(:,:,2:end).*coefs);
   Q.d = P.d+n;
end;

props(Q,P.v,P.h);

%end .. @pol/deriv
