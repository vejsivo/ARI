function Pshift = shift(P,n,symb);
%SHIFT  Shift polynomial
%
% PS = SHIFT(P,N) with polynomial  matrix P and integer scalar N
% computes PS(VAR) = P(VAR)*VAR^N, where VAR is the variable of P.
% If P has no variable and the variable for the result is needed,
% the global variable is used.
%
% PS = SHIFT(P,N,SYMB) act as above but the variable symbol of 
% the result is set to SYMB.
%
% N can be also an integer matrix of the same dimensions as P.
% The shift is performed elementwise: every P(i,j) is shifted
% by N(i,j). In the case when N is matrix, P can be also scalar;
% it is shifted by every entry of matrix N, the dimensions of
% the result being the same as those of N.
%
% With some entry of N negative, the result may contain negative
% powers of VAR. In such a case, the result is a two-sided
% polynomial, its symbol being always 'z'. When VAR is neither
% 'z' nor 'z^-1', a warning is issued and VAR is changed before
% the shift is made: 'q','s' or 'p' to 'z' but 'd' to 'z^-1'.
%
% See also POL/TIMES, POL/MTIMES, POL/POWER, POL/MPOWER.

%        Author(s): M. Hromcik, M. Sebek 20-5-98
%        Copyright (c) 1998 by Polyx, Ltd.
%        $Revision: 2.0 $  $Date: 10-Jun-1998  10:28:34  $
%        $Revision: 3.0 $  $Date: 06-Apr-2000  J. Jezek  $
%                          $Date: 06-Feb-2001  J. Jezek  $
%                          $Date: 09-Jul-2001  J. Jezek  $ 

global PGLOBAL;

na = nargin;
if na>=2,
   if ~isa(n,'double') | ndims(n)>2 | ~isreal(n) | ...
         (~isempty(n) & any(any(floor(n)~=n)) ),
      error('Invalid shift; must be integer.');
   end;
else
   error('Not enough input arguments.');
end;

if na==3 & ~isempty(symb),
   if isa(symb,'pol'),
         [vs1,vs2,vd] = size(symb);
         if all([vs1,vs2,vd]==1)&(~any(symb.c(:,:)-[0,1])),
            symbh = symb.h; symb = symb.v;
         else
            error('Invalid variable symbol.');
         end;
   elseif ~ischar(symb),
      error('Invalid variable symbol.');
   end;
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
   symb = P.v; symbh = P.h;   
end; 
if isempty(symb),
   symb = PGLOBAL.VARIABLE;
   var = pol([0 1],1,symb); symbh = var.h;
end;

eval('P = pol(P);', 'error(peel(lasterr));');
ns = size(n); Ps = P.s;
if all(ns==1),           % n  is scalar
   if any(Ps~=1),                  % P  is matrix
      n = ones(Ps)*n;
   end;
else                     % n  is matrix
   if all(Ps==1),                  % P  is scalar
      P = ones(ns)*P; Ps = P.s;
   else                            % P  is matrix
      if any(Ps~=ns),
         error('Matrices not of the same dimensions.');
      end;
   end;
end;

if isempty(P.c),
   Pshift = P; return;
end;

nmin = min(min(n));
if nmin<0,
   for i = 1:Ps(1),
      for j = 1:Ps(2),
         nij = n(i,j);
         if nij<0,
            if ~isempty(nonzeros(P.c(i,j,1:min(abs(nij),end)))),
               if strcmp(P.v,'z^-1') | ...
                     (isempty(P.v) & strcmp(symb,'z^-1')),
                  n = -n;
               elseif strcmp(P.v,'d') | ...
                     (isempty(P.v) & strcmp(symb,'d')),
                  n = -n; P.v = 'z^-1';
                  warning('Variable symbol changed to ''z^-1''.');
               elseif ~isempty(P.v) & ~strcmp(P.v,'z'),
                  P.v = 'z';
                  warning('Variable symbol changed to ''z''.');
                  if symbh==0, symbh = NaN;
                  end;
               end;
               P.h = symbh;
               Pshift = shift(tsp(P),n);
               return;
            end;
         end;
      end;
   end;
end;

nmax = max(max(n));
Psd = P.d + nmax;
Pshift.d = Psd;
Pshift.s = Ps;
Pshift.c = zeros(Ps(1),Ps(2),Psd+1);
Pshift.v = symb;
Pshift.h = symbh;
Pshift.u = [];
Pshift.version = 3.0;
Pshift = class(Pshift,'pol');

Pd1 = P.d + 1;
for i = 1:Ps(1),
   for j = 1:Ps(2),
      nij = n(i,j);
      if nij>=0,
         Pshift.c(i,j,1+nij:Pd1+nij) = P.c(i,j,1:Pd1);
      else
         Pshift.c(i,j,1:Pd1+nij) = P.c(i,j,-nij+1:Pd1);
      end;
   end;
end;
Pshift = pclear(Pshift);

%end .. @pol/shift
