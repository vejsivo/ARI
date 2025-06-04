function Y = mvalue(F,X)
%MVALUE      Matrix value of two-sided polynomial
%
% The command  Y = MVALUE(F,X)  returns value F(X)
% of two-sided polynomial function of the argument X.
%
% If F is scalar tsp then X may be square matrix;
% otherwise it must be scalar.
%
% When F is empty or constant, the function does not
% depend on X. In such a case, X may be any scalar or
% any compatible matrix or it may be omitted.
%
% If X is scalar it may be zero or infinite.
%
% See also POL/MVALUE, TSP/VALUE.

%     Author: J.Jezek, 09-Oct-2000
%     Copyright(c) 2000 by Polyx, Ltd.
%     $ Revision $  $ Date 31-Oct-2000 $
%                   $ Date 19-Dec-2000 $
%                   $ Date 27-Apr-2003 $

ni = nargin;
if ni==2,
   if ~isa(X,'double') | ndims(X)~=2,
      error('Invalid 2nd argument.');
   end;
end;

if ni==1,
   [FD,FC] = declass(F);
   if strcmp(FC,'double'),
      Y = FD; return;
   else
      error('Not enough input arguments.');
   end;
end;

if length(X)==1,
   if isinf(X),
      Y = mvalue(nneg(F),X); return;
   elseif X==0,
      Y = mvalue(npos(F),Inf); return;
   end;
end;

Y = 0;
eval('Y = mvalue(F.p,X);',...
   'error(peel(lasterr));');
Y = X^F.o * Y;

%end .. @tsp/mvalue
