function P = pzer(P, zeroing)
%PZER  Perform zeroing on a two-sided polynomial matrix P
%
% P = PZER(P) sets the coefficients in P whose absolute value is less
% than the current value of global zeroing tolerance equal to zero.
%
% P = PZER(P,ZEROING) uses the value of argument ZEROING as the 
% local tolerance.
% 
% To suppres zeroing globally set the global zeroing tolerance equal to zero.

%       Author(s):  S. Pejchova  24-7-00
%       Copyright (c) 2000 by Polyx, Ltd.
%       $Revision: 3.0 $  $Date: 25-Jul-2000 16:28:34  S.Pejchova  $


global PGLOBAL;
eval('in_tol=PGLOBAL.ZEROING;',...
      'error(''Use PINIT to initialize Polynomial Toolbox.'');');
if nargin==2,
   in_tol=zeroing;
elseif nargin < 1,
   error('Not enough input arguments.');
end;

if in_tol > 0, % if ZEROING <= 0, zeroing is not performed
   P.p=pzer(P.p,in_tol);
end;

% second, remove zero leading coefficients
% that may have been introduced
P = tclear(P);

%end .. @tsp/pzer
