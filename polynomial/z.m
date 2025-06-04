function z = z(T);
%Z  Create a simple basic polynomial (monomial)
%
% Z is used to enter the first power of the variable z.
% It is interpreted as the discete time shift (advance) operator.
%
% Since Z is a function it may be overridden and used as a variable.
%
% Optional argument T may specify sampling period. Default is 1.
%
% See also S, P, Q, D, ZI, V.

%       Author(s):  S. Pejchova, M. Sebek 17-4-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $ Revision: 2.0 $  $ Date: 25-Jun-1998 12:35:34   $
%       $           3.0 $  $ Date: 16-Dec-2000 J. Jezek   $
%                          $ Date: 15-Jul-2001 J. Jezek   $
%
% Effect on other properties:
% Variable is set equal to z.
% UserData are deleted.

global PGLOBAL;
eval('PGLOBAL.VARIABLE;', 'painit;');

z = pol([0 1],1,'z');
if nargin==1,
   if isa(T,'double') & (isempty(T) | ...
         (length(T)==1 & isreal(T) & T>=0)),
      eval('z.h = T;', 'error(peel(lasterr));');
   else
      error('Invalid sampling period.');
   end;
end;

%end .. z
