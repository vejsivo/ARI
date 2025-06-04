function [b,blog] = integral(a,var)
%INTEGRAL      Integral of polynomial
%
% The command  B = INTEGRAL(A,VAR)  returns the integral
% of polynomial A. The integration constant is selected
% such that the absolute term of B is zero.
%
% The optional input argument VAR is the integration
% variable. For polynomial in s,p,q,d, this argument,
% if present, must be the variable symbol of A.
% For polynomials in z or z^-1, it may be any of
% these two; then the integral may happen to be
% a two-sided polynomial or to contain a logarithmic
% term. In the latter case, this term can be obtained
% as the logarithm of the second output argument.
%
% See also POL/DERIV.

%     Author: J.Jezek, 05-Oct-2000
%     Copyright(c) by Polyx, Ltd.
%     $ Revision $  $ Date 01-Aug-2001  sampl per $
%                   $ Date 22-Sep-2002            $

global PGLOBAL;

ni = nargin;
if ni<1,
   error('Not enough input arguments.');
end;
eval('a = pol(a);', ...
   'error(''Invalid 1st argument.'');');
s1 = a.s(1); s2 = a.s(2);
ons = ones(s1,s2); zers = zeros(s1,s2);
no = nargout;
if no==2,
   blog = pol(zers);
end;
h = [];

if ni==2,
   if isa(var,'char'),
      if isempty(var),
         var = PGLOBAL.VARIABLE;
      elseif isempty(strmatch(var,{'z';'z^-1';'zi';'s';'p';'d';'q'},'exact')),
         error('Invalid integration variable.');
      end;
      w = pol([0 1],1,var);
   else
      eval('w = pol(var);', ...
         'error(''Invalid integration variable.'');');
     [vs1,vs2,vd] = size(w);
      if all([vs1,vs2,vd]==1) & all(var.c(:,:)==[0,1]),
         h = w.h; var = w.v;
      else
         error('Invalid integration variable.');
      end;
   end;
else
   var = PGLOBAL.VARIABLE;
   w = pol([0 1],1,var);
end;

if isempty(a.v),
   b = a*w;
   
else
   direct = 1;
   if ni==2,
      if strcmp(var,'zi'), var = 'z^-1';
      end;
      if ~strcmp(var,a.v),
         if (strcmp(a.v,'z') & strcmp(var,'z^-1')) | ...
               (strcmp(a.v,'z^-1') & strcmp(var,'z')),
            direct = 0;
         else
            error('Invalid integration variable.');
         end;
      end;
   else
      var = a.v;
   end;
   
   if ~isempty(h) & isfinite(h) & ~isempty(a.h) & isfinite(a.h),
      if a.h~=h,
         warning('Inconsistent sampling periods.');
         a.h = NaN;
      end;
   end;
   
   if isempty(a.h) & ~isempty(h),
      a.h = h;
   end;
   
   d = a.d;
   if direct,      % integration variable == a.v
      dp1 = d+1;
      invnum = 1./(1:dp1);
      kronnum = kron(invnum,ons);
      b.d = dp1;
      b.s = a.s;
      b.c = zeros(a.s(1),a.s(2),dp1+1);
      b.c(:,:) = [zers, a.c(:,:) .* kronnum];
      
   else            % integration variable == a.v^-1
      w = pol([0 1],1,var);
      w.h = a.h;
      ac2 = a.c(:,:,2);
      if (any(any(ac2))),
         if no<=1,
            warning('Integral is not polynomial, contains logarithm.');
         else
            blog = ac2*w;
         end;
      end;
      dm1 = d-1;
      invnum = 1./(1:dm1);
      kronnum = kron(invnum,ons);
      b.d = dm1;
      b.s = a.s;
      b.c = zeros(a.s(1),a.s(2),d);
      ac3 = a.c(:,:,3:d+1);
      b.c(:,:) = [zers, -ac3(:,:).*kronnum];
   end;
   
   b.v = a.v;
   b.h = a.h;
   b.u = [];
   b.version = 3.0;
   b = class(b,'pol');
   b = pclear(b);
   
   if ~direct,
      ac1 = a.c(:,:,1);
      if (any(any(ac1))),
         b = b + ac1*w;
      end;
   end;
   
end;

%end .. @pol/integral

      
      