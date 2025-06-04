function P = pclear(P)
%PCLEAR  Perform clearing on a polynomial matrix.
% 
%   This macro will remove extra zero leading coefficients in
%   polynomial matrix P, if any.

%   This macro is most notably called in macro POL/SUBSREF when
%   extracting a submatrix, and in macro PZER when performing zeroing.

%       Author(s):  M. Hromcik, S. Pejchova, M. Sebek 12-2-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 3.0 $  $Date: 08-Nov-1999 14:51:34   $
%                         $Date: 18-May-2000  J. Jezek  $
%                         $Date: 26-Jul-2000  J. Jezek, S. Pejchova  $
%                         $Date: 01-Aug-2001  J. Jezek  $
%                         $Date: 15-Sep-2009  M. Sebek  $ modified NaN
%                         handling


if ~isempty(P.c),
%   P.c = P.c(:,:,1:max(find(any(any(P.c,1),2)))); original version
    P.c =  P.c(:,:,1:max(find(any(any( P.c,1),2)+any(any(isnan( P.c),1),2))));
   P.d = size(P.c,3) - 1;
else
   P.v = '';
end;
if P.d <= 0,
   P.v = '';
   if P.d < 0, P.d = -Inf;
   end;
end;
if isempty(P.v), P.h = [];
end;

%end .. @pol/private/pclear
