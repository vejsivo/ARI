function B = dilate(A,k,h,var)
%DILATE     Dilate polynomial
%
% Let A be a polynomial in variable 'z^-1' or 'd'
%    A = A(0) + ... + A(n)*z^-n + ...
% whose n-th coefficient is A(n). The command
%    B = DILATE(A,K,H)
% where scalar integers K,H are dilating period, K>=1, and dilating
% phase, H>=0, returns polynomial B(z^-1) = z^-H * A(z^-K). Note
% that only such coefficients of B are nonzero which can be expressed
% B(K*n+H) = A(n) with some n. The arguments K,H are optional,
% the defaults being  K=2, H=0.
%
% For a polynomial in variable 'z' or 'q', the numbering of H is
% reversed, the formulas being  B(z) = z^-H * A(z^K) ,  B(K*n-H) = A(n). 
% This rule is for compatibility with TSP/DILATE. For H~=0, the result
% may be a two-sided polynomial or a fraction.
%
% For a polynomial in variable 's' or 'p', the dilatation is
% not defined.
%
% The variable symbol of B is taken from that of A. If A is constant
% and has no symbol then the symbol of B may be taken from an optional
% (fourth) argument VAR or from the global discrete time variable symbol.
%
% The meaning of the variable symbol has changed:
%  old var = new var to the K-th power.
% The sampling period of B is K-times less than that of A. When VAR
% is applied, the sampling period of B is K-times less than that of VAR.
%
% See also TSP/DILATE, POL/RESAMP.

%        Author:  J. Jezek 08-Nov-1999
%        Copyright(c) 1999 by Polyx, Ltd.
%        $ Revision $  $ Date 02-Jun-2000 $
%                      $ Date 04-Oct-2000 $
%                      $ Date 30-Sep-2002 $
%                      $ Date 28-Feb-2003 $

global PGLOBAL;

ni = nargin;
if ni<1,
   error('Not enough input arguments.');
end;
eval('A = pol(A);', ...
   'error(''Invalid 1st argument.'');');

if ni<2 | isempty(k),
   k = 2;
else
   if ~isa(k,'double') | length(k)~=1 | ~isreal(k) | ...
         floor(k)~=k | k<1,
      error('Invalid dilating period.');
   end;
end;

if ni<3 | isempty(h),
   h = 0;
else
   if ~isa(h,'double') | length(h)~=1 | ~isreal(h) | floor(h)~=h ...
         | h<0,
      error('Invalid dilating phase.');
   end;
end;

Av = A.v;
if strcmp(Av,'s') | strcmp(Av,'p'),
   error('Invalid variable symbol; must be discete time.');
elseif strcmp(Av,'z') | strcmp(Av,'q'),
   h = -h;
end;

if isempty(A) | all(all(A==0)),
   B = A; return;
end;

if h>=0,
   B.d = A.d*k + h;
   B.s = A.s;
   if ~isempty(A.c),
      B.c = zeros(A.s(1),A.s(2),B.d+1);
      B.c(:,:,1+h:k:end) = A.c;
   else
      B.c = zeros(A.s(1),A.s(2),0);
   end;
   per = []; varused = logical(0);
   if ~isempty(Av),
      B.v = Av; per = A.h;
   elseif B.d>0,
      if ni==4 & ~isempty(var),
         if isa(var,'char'),
            per = 1;
         else
            eval('pvar = pol(var);', ...
               'error(''Invalid 4th argument.'');');
            [vs1,vs2,vd] = size(pvar);
            if all ([vs1,vs2,vd]==1) & all(pvar.c(:,:)==[0 1]),
               var = pvar.v; per = pvar.h;      
            else error('Invalid 4th argument.');
            end;
         end;
         if isempty(strmatch(var,{'z';'z^-1';'zi';'q';'d'},'exact')),
            error('Invalid 4th argument.');
         end;
      else
         var = PGLOBAL.DISCRVAR; per = 1;
      end;
      B.v = var; varused = logical(1);
   else
      B.v = '';
   end;
   B.h = per;
   if ~isempty(B.h), B.h = B.h/k;
   end;
   B.u = [];
   B.version = 3.0;
   B = class(B,'pol');
   if varused & (strcmp(B.v,'z') | strcmp(B.v,'q')),
      Bv = B.v; B.v = 'z^-1';
      B = sdf(B); B = reverse(B);
      B.v = Bv;
   end;
else
   B.d = A.d*k;
   B.s = A.s;
   B.c = zeros(A.s(1),A.s(2),B.d+1);
   if ~isempty(A.c),
      B.c(:,:,1:k:end) = A.c;
   end;
   B.v = Av;
   B.h = A.h;
   if ~isempty(B.h), B.h = B.h/k;
   end;
   B.u = [];
   B.version = 3.0;
   B = class(B,'pol');
   if strcmp(Av,'z'),
      B = shift(B,h);
   else   %strcmp(Av,'q')
      qh = q^h; qh.h = B.h;
      B = qh*B;
   end;
end;

%end .. @pol/dilate
   