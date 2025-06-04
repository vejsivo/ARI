function q = q(T);
%Q  Create a simple basic polynomial (monomial)
%
% Q is used to enter the first power of the variable q.
% It is interpreted as the discrete time shift (advance) operator.
%
% Since Q is a function it may be overridden and used as a variable.
%
% Optional argument T may specify sampling period. Default is 1.
%
% See also S, P, Z, D, ZI, V.

%       Author(s):  S. Pejchova, M. Sebek 17-4-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 25-Jun-1998 12:31:34   $
%       $          3.0 $  $Date: 16-Dec-2000 J.Jezek    $
%                         $Date; 15-Jul-2001 J.Jezek    $
%
% Effect on other properties:
% Variable is set equal to q.
% UserData are deleted.

global PGLOBAL;
eval('PGLOBAL.VARIABLE;', 'painit;');

q = pol([0 1],1,'q');
if nargin==1,
   if isa(T,'double') & (isempty(T) | ...
         (length(T)==1 & isreal(T) & T>=0)),
      eval('q.h = T;', 'error(peel(lasterr));');
   else
      error('Invalid sampling period.');
   end;
end;

%end .. q
