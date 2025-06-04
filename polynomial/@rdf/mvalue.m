function Y = mvalue(F,X)
%MVALUE      Matrix value of right-denominator fraction
%
% The command  Y = MVALUE(F,X)  returns value F(X)
% of right-den fractional function of the argument X.
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
% See also POL/MVALUE, POL/VALUE, RDF/VALUE.

%     Author: J.Jezek, 09-Oct-2000
%     Copyright(c) 2000 by Polyx, Ltd.
%     $ Revision $  $ Date 01-Nov-2000 $
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
   dF = max(deg(num),deg(den));
   N = num{dF}; D = den{dF};
else
   F = coprime(F); N = 0;
   eval('N = mvalue(F.frac.num,X);',...
      'error(peel(lasterr));');
   D = mvalue(F.frac.den,X);
end;
Y = N/D;

%end .. @rdf/mvalue
