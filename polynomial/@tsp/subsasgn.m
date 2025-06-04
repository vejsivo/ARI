function T = subsasgn(T,S,R)
%SUBSASGN   Subscripted assignment for two-sided polynomial
%
% If T is the two-sided polynomial matrix 
%   TNK*z^-k  + ... + TN1*z^-1 + T0 + TP1*z + TP2*z^2 + .. + TPM*z^m
%
% then
%    T{-i:j}  = R  is  [TNI ... TN1 T0 TP1 TP2 .. TPJ] = R,
%    T{k}     = R  is the k-term of T, TK = R.
%
%    T(ri,ci) = R  sets the polynomial submatrix of T of row indices ri
%                  and column indices ci equal to R. 
%
%    T.h      = R  sets the sampling period of T.
%
%    T.user   = R  sets T.u, the user data of T.
%
% The above subscripted references may be combined, i.e.,
% T{di}(ri, ci) sets the constant coefficient matrix
% corresponding to a submatrix T equal to R.
%
% See also TSP/SUBSREF.

%       Author(s):  S. Pejchova  28-07-99
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 3.0 $  $Date: 27-Mar-2000 14:48:34   $
%                         $Date: 22-May-2000 10:00:00  J.Jezek  $
%                         $Date: 28-Feb-2002  J.Jezek  $

if nargin<3,
   error('Not enough input arguments.');
end;
if ~isstruct(S),
   error('Invalid 2nd argument.');
end;
eval('T = tsp(T);', 'error(''Invalid 1st argument.'');');

S_lg = length(S);
if S_lg==0, return;
end;
St1 = S(1).type;
switch St1,
case '.',
   Ss1 = S(1).subs; flstr = lower(Ss1(1));
   switch flstr,
   case 'v',
      if S_lg>1,
         error('Invalid combination in subscripted asignment.');
         end;
      ok = 0;
      if isa(R,'pol'),
         [rs1,rs2,rd] = size(R);
         if all([rs1,rs2,rd]==1 & ~any(R.c(:,:)-[0 1])),
            R = R.v;
         end;
      end;
      if isa(R,'char'),
         if ~isempty(R),
            if strcmp(R,'z') | strcmp(R,'zi') | strcmp(R,'z^-1'),
               ok = 1;
            end;
         else
            Td = T.d;
            if isempty(Td) | isinf(Td) | T.t==Td,
               ok = 1;
            end;
         end;
      end;
      if ~ok, error('Invalid variable symbol.');
      end;
      
   case 'h',
      if S_lg>1,
         eror('Invalid combination in subscripted asignment.');
      end;
      if isempty(R) | isempty(T.v),
         T.h = [];
      elseif ~isa(R,'double') | length(R)~=1 | ...
            ~isreal(R) | R<0 | isinf(R),
         error('Invalid sampling period.');
      else T.h = R;
      end;
      eval('T.p.h = R;', 'error(peel(lasterr));');
      T = tclear(T);
      
   case 'u',
      S(1) = [];
      if isempty(S), T.u = R;
      else
         Tu = T.u;
         eval('Tu = subsasgn(Tu,S,R);', ...
            ['if isa(T.u,''pol'') | isa(T.u,''tsp'') | isa(T.''frac''),', ...
               'error(peel(lasterr)); else error(peelf(lasterr)); end;']);
         T.u = Tu;
      end;
      
   otherwise
      error('Invalid field in subscripted reference.');
   end;
   
case '()',
   eval('R = tsp(R);', ...
      'error(''Argument canot be assigned to tsp.'');');
   [th,Th,T,R] = testht(T,R);
   if th==0,
      warning('Inconsistent sampling periods.');
   end;
   Tu = T.u;
   oo = min(T.o,R.o);
   TT = shift(T.p,T.o-oo,'z');
   RR = shift(R.p,R.o-oo,'z');
   eval('TT = subsasgn(TT,S,RR);', ...
      'error(peel(lasterr));');
   T = tsp(TT); T.h = Th; T.u = Tu;
   T = shift(T,oo);

case '{}',
   Th = T.h; Tu = T.u;
   Ss1 = S(1).subs;
   if length(Ss1)>1,
      error('More than one subscript in {} .');
   end;
   k = Ss1{1};
   if isa(k,'double'),
      oo = min(min(T.o,k));
      S(1).subs{1} = k-oo;
   else
      oo = T.o;
   end;
   TT = shift(T.p,T.o-oo,'z');
   eval('TT = subsasgn(TT,S,R);','error(peel(lasterr));');
   TT.v = 'z'; T = tsp(TT);
   T.h = Th; T.u = Tu;
   T = shift(T,oo);
      
end;

%end .. @tsp/subsasgn

      