function v = v(T);
%V  Create a simple basic polynomial (monomial) 
%   with the current value of the global variable string
%
% V is used to enter the first power of the current value 
% of the global variable string.
%
% Since V is a function it may be overridden and used as a variable.
%
% For a discrete time variable, a smpling period may be
% specified as an optional argument. The default is 1.
% For a continuous time variable, the sampling period is always 0.
%
% See also S, P, Z, Q, D, ZI.

%       Author(s):  S. Pejchova, M. Sebek 17-4-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 25-Jun-1998 12:34:34   $
%       $Revision: 3.0 $  $Date: 18-Jul-2000 10:19:34  S.Pejchova  $
%                         $Date: 15-Jul=2001 J.Jezek  $
%
% Effect on other properties:
% Variable is set equal to global variable string.
% UserData are deleted.

global PGLOBAL;
eval('PGLOBAL.VARIABLE;', 'painit;');

v = pol([0 1],1);
if nargin==1,
   if isa(T,'double') & (isempty(T) | ...
         (length(T)==1 & isreal(T) & T>=0)),
      eval('v.h = T;', 'error(peel(lasterr));');
   else
      error('Invalid sampling period.');
   end;
end;

%end .. v
