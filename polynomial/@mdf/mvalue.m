function Y = mvalue(F,X)
%MVALUE      Matrix value of matrix-denominator fraction
%
% The command  Y = MVALUE(F,X)  returns value F(X)
% of matrix-den fractional function of the argument X.
%
% If F is scalar fraction then X may be square
% matrix; otherwise it must be scalar.
%
% When F is empty or constant, the function does not
% depend on X. In such a case, X may be any scalar or
% any compatible matrix or it may be omitted.
%
% If X is scalar, it may be infinite.
%
% See also POL/MVALUE, POL/VALUE, MDF/VALUE.

%     Author: J.Jezek, 09-Oct-2000
%     Copyright(c) 2000 by Polyx, Ltd.
%     $ Revision $  $ Date 06-Nov-2000 $
%                   $ Date 16-Dec-2000 $
%                   $ Date 14-Oct-2002 $
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

if isempty(F),
   Y = zeros(size(F)); return;
end;

if length(X)==1 & isinf(X),
   num = F.frac.num; den = F.frac.den;
   dF = max(deg(num,'ent'),deg(den,'ent'));
   sF = size(F); sF1 = sF(1); sF2 = sF(2);
   N = zeros(sF); D = N;
   for i = 1:sF1,
      for j = 1:sF2,
         dFij = dF(i,j);
         N(i,j) = num{dFij}(i,j);
         D(i,j) = den{dFij}(i,j);
      end;
   end;
else
   N = 0;
   eval('N = mvalue(F.frac.num,X);',...
      'error(peel(lasterr));');
   D = mvalue(F.frac.den,X);
end;

if length(F)==1,
   Y = N/D;
else
   Y = N./D;
end;

%end .. @mdf/mvalue
