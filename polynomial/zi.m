function zi = zi(T);
%ZI  Create a simple basic polynomial (monomial)
%
% ZI is used to enter the first power of the variable z^-1.
% It is interpreted as the discrete time shift (delay) operator.
%
% Since ZI is a function it may be overridden and used as a variable.
%
% Optional argument T may specify sampling period. Default is 1.
%
% See also S, P, Z, Q, D, V.

%       Author(s):  S. Pejchova, M. Sebek 17-4-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 25-Jun-1998 12:38:34   $
%       $          3.0 $  $Date: 16-Dec-2000 J.Jezek    $
%                         $Date: 15-Jul-2001 J.Jezek    $
%
% Effect on other properties:
% Variable is set equal to z^-1.
% UserData are deleted.

global PGLOBAL;
eval('PGLOBAL.VARIABLE;', 'painit;');

zi = pol([0 1],1,'z^-1');
if nargin==1,
   if isa(T,'double') & (isempty(T) | ...
         (length(T)==1 & isreal(T) & T>=0)),
      eval('zi.h = T;', 'error(peel(lasterr));');
   else
      error('Invalid sampling period.');
   end;
end;

%end .. zi
