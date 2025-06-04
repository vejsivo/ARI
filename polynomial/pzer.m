function Pz = pzer(P, zeroing)
%PZER  Perform zeroing on a polynomial matrix
%
% PZ = PZER(P) sets those coefficients in P equal to zero
% whose absolute values are less than the current value of 
% the global zeroing tolerance.
%
% PZ = PZER(P,ZEROING) uses the value of input argument 
% ZEROING as the local tolerance.
% 
% To suppress zeroing globally set the global zeroing 
% tolerance equal to zero.

%       Author(s):  D. Henrion, S. Pejchova, M. Sebek 24-2-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 25-Jun-1998 11:45:34   $
%       $Revision: 3.0 $  $Date: 25-Jul-2000 16:26:34  S.Pejchova  $
%                         $Date: 13-Jul-2001  J.Jezek  $
%  No effect on other properties.

if nargin < 1,
   error('Not enough input arguments.');
elseif nargin==1,
   eval('Pz=pzer(pol(P));', 'error(peel(lasterr));');
else,
   eval('Pz=pzer(pol(P),zeroing);', 'error(peel(lasterr));');
end;

%end .. pzer
