function p = p(T);
%P  Create a simple basic polynomial (monomial)
%
% P is used to enter the first power of the variable s.
% It is interpreted as the continuous time derivative operator.
%
% Since P is a function it may be overridden and used as a variable.
%
% See also S, Z, Q, D, ZI, V.

%       Author(s):  S. Pejchova, M. Sebek 17-4-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 25-Jun-1998 12:24:34   $
%       $Revision: 3.0 $  $Date: 18-Jul-2000 10:11:34  S.Pejchova  $
%                         $Date; 15-Jul-2001 J.Jezek $
%
% Effect on other properties:
% Variable is set equal to p.
% UserData are deleted.

global PGLOBAL;
eval('PGLOBAL.VARIABLE;', 'painit;');

p = pol([0 1],1,'p');
if nargin==1,
   if isa(T,'double') & (isempty(T) | ...
         (length(T)==1 & isreal(T) & T>=0)),
      eval('p.h = T;', 'error(peel(lasterr));');
   else
      error('Invalid sampling period.');
   end;
end;

%end .. p
