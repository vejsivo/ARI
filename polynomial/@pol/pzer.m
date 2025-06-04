function P = pzer(P, zeroing)
%PZER  Perform zeroing on a polynomial matrix P
%
% P = PZER(P) sets the coefficients in P whose absolute value is less
% than the current value of global zeroing tolerance equal to zero.
%
% P = PZER(P,ZEROING) uses the value of argument ZEROING as the 
% local tolerance.
% 
% To suppres zeroing globally set the global zeroing tolerance equal to zero.

%       Author(s):  D. Henrion, S. Pejchova, M. Sebek 24-2-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 25-Jun-1998 11:34:34   $
%       $Revision: 3.0 $  $Date: 21-Jul-2000 12:20:34  S.Pejchova  $

%  No effect on other properties.

global PGLOBAL;
eval('in_tol=PGLOBAL.ZEROING;', 'painit;');
if nargin == 2,
   if ~isa(zeroing,'double') | length(zeroing)~=1 | ...
         ~isreal(zeroing),
      error('Invalid tolerance.');
   end;
   in_tol=zeroing;
elseif nargin < 1,
   error('Not enough input arguments.');
end;

if in_tol > 0, % if ZEROING <= 0, zeroing is not performed
   % first, set neglected entries to zero
   if isreal(P.c),
      P.c(abs(P.c)<in_tol)=0;
   else,
      Pr=real(P.c); Pi=imag(P.c);
      Pr(abs(Pr)<in_tol)=0;
      Pi(abs(Pi)<in_tol)=0;
      P.c=Pr+i*Pi;
   end;
end;

% second, remove zero leading coefficients
% that may have been introduced
P = pclear(P);

%end .. @pol/pzer
