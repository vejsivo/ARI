function T = subsref(T,S)
%SUBSREF  Subscripted reference for two-sided polynomial
%
% If T is the polynomial matrix 
%   TNK*z^-k  + ... + TN1*z^-1 + T0 + TP1*z + TP2*z^2 + .. + TPM*z^m
%
% then
%    T{di} 
% returns the block row constant matrix corresponding to 
% the degree indices di, in particular
%    T{-i:j} is [TNI ... TN1 T0 TP1 TP2 .. TPJ],
%    T{0} is T0, the constant term of T, and
%    T{k} is Tk, the  k-term of T.
%    T(ri, ci) returns the polynomial submatrix of T of row 
%              indices ri and column indices ci.
%    T(X) where X(var) is scalar polynomial, tsp or fraction,
%         returns the composition  T(X(var)).
%         If possible, the result is tsp;
%         if not, it is scalar-denominator fraction. 
%
%    T.pol returns the POL object  corresponding to the TSP object,
%          but shifted to positive powers of 'z'.
%    T.size returns the size of T.
%    T.deg returns the degree of T.
%    T.tdeg returns the trailing degree of T.
%    T.offset is a offset of T.pol (the degree of 'z' to obtain 
%          from T.pol the original)
%    T.range returns the degree of the corresponding POL object T.p.
%    T.var returns the indeterminate string variable of T.
%    T.h  returns the sampling period.
%    T.user returns the user data of T.
%
% The above subscripted references may be combined, i.e.
% T(ri, ci).deg returns the degree of a submatrix of T, and
% T{di}(ri, ci) returns the constant coefficient matrix
% corresponding to the degree di of the submatrix T(ri,ci).
%
% See also TSP/SUBSASGN.

%       Author(s):  S. Pejchova  21-07-99
%       Copyright (c) 1999 by Polyx, Ltd.
%       $Revision: 3.0 $  $Date: 17-Sep-1999 15:01:34  $
%                         $Date: 22-May-2000  J.Jezek  $
%                         $Date: 31-Oct-2000  J.Jezek  $
%                         $Date: 01-Feb-2002  J.Jezek  $
%                         $Date: 03-Jul-2002  J.Jezek  $ 
if nargin<2,
   error('Not enough input arguments.');
end;
if ~isstruct(S),
   error('Invalid 2nd argument.');
end;

To = T.o; Tp = T.p;
S_lg = length(S); S1 = S(1);
St1 = S1.type; Ss1 = S1.subs;
Ss1_lg = length(Ss1);
switch St1,
case '.',
   flstr = lower(Ss1); done = 0;
   if length(flstr)>=3,
      flstr = flstr(1:3);
      if strcmp(flstr,'ver'),
         T = T.version; done = 1;
      elseif strcmp(flstr,'den'),   % denominator
         if To>=0, T = pol(1);
         else T = z^-To;
         end;
         done = 1;
      end;
   end;
   if ~done,
      flstr1 = flstr(1);
      switch flstr1,
      case 'f',
         if S_lg==1,
            error('Invalid field in subscripted reference.');
         end;
         S(1) = []; S_lg = length(S);
         St1 = S(1).type; Ss1 = S(1).subs;
         if ~strcmp(St1,'.'),
            error('Invalid combination in subscripted reference.');
         end;
         flstr = lower(Ss1); flstr1 = flstr(1);
         switch flstr1,
         case 'n',
            if To>=0, T = pol(T);
            else T = Tp;
            end;
         case 'd',
            if To>=0, T = pol(1);
            else T = z^-To;
            end;
         otherwise
            error('Invalid field in subscripted reference.');
         end;
      case 'n',
         if To>0, T = pol(T);
         else T = Tp;
         end;
      case 'p', T = T.p;
      case 's', T = T.s;
      case 'd', T = T.d;
      case 't', T = T.t;
      case 'o', T = T.o;
      case 'r', T = T.r;
      case 'v', T = T.v;
      case 'h', T = T.h;
      case 'u', T = T.u;
      otherwise,
         error('Invalid field in subscripted reference.');
      end;
   end;
   
case '()',
   if Ss1_lg==1,
      X = Ss1{1};
      if isa(X,'pol') | isa(X,'tsp') | isa(X,'frac'),
         Y = 0;
         eval('Y = subsref(sdf(T),S);','error(peel(lasterr));');
         if isa(Y,'frac') & (strcmp(Y.v,'z') | strcmp(Y.v,'z^-1')),
            eval('T = tsp(Y);','T = Y;');
         else
            T = Y;
         end;
         return;
      else
         eval('PP = Tp(X);','error(peel(lasterr));');
         T = tsp(PP); T = shift(T,To);
      end;
   elseif Ss1_lg==2,
      X = Ss1{1}; Y = Ss1{2};
      eval('PP = Tp(X,Y);','error(peel(lasterr));');
      T = tsp(PP); T = shift(T,To);
   else
      error('More than two subscripts.');
   end;
   
case '{}',
   if Ss1_lg==1,
      X = Ss1{1};
      if ~isfinite(X),
         error('Subscript is not finite.');
      end;
      if X<To, T = zeros(size(T));
      else
         eval('T = Tp{X-To};','error(peel(lasterr));');
      end;
   else
      error('More than one subscript in {}.');
   end;
   
otherwise,
   error('Invalid 2nd argument.');
end;

S(1) = [];
if ~isempty(S),
   if (isa(T,'double') | isa(T,'char')) & strcmp(S(1).type,'.'),
      error('Invalid combination in subscripted reference.');
   end;
   eval ('T = subsref(T,S);', ...
      ['if isa(T,''pol'') | isa(T,''tsp'') | isa(T,''frac''),', ...
         'error(peel(lasterr)); else error(peelf(lasterr)); end;']);
end;

%end .. @tsp/subsref

         