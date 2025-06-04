function T = tclear(T)
%TCLEAR  Perform clearing on a two-sided polynomial matrix.
% 
%   This macro will remove extra zero leading and trailing
%   coefficients, it means the zero matrix coefficients at
%   lowest and highest degrees of 'z' of TSP object T,
%   if any.

%       Author(s):  S. Pejchova  13-07-99
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 3.0 $  $Date: 17-Sep-1999 11:16:34   $
%                         $Date: 13-Oct-1999 12:00:00  J.Jezek  $
%                         $Date: 22-May-2000 16:00:00  J.Jezek  $
%       $Revision: 3.0 $  $Date: 24-Jul-2000 14:01:34  S.Pejchova  $

Tr = T.p.d; Tc = T.p.c;
if ~isempty(Tc),
   offst=T.o; 
   fT = find(any(any(Tc,1),2));
   Tc = Tc(:,:,min(fT):max(fT));
   Tc3 = size(Tc,3);
   T.r = Tc3-1;
   T.o = offst+(Tr+1-Tc3);
   T.d = offst+Tr;
   T.t = T.o;
   T.p = pol(Tc(:,:),Tc3-1,'z');
   T.p.h = T.h;
   Tr = T.r;
elseif Tr < 0,
   T.o=0; T.t=inf; T.d=-inf; T.r=-inf;
end;
if isempty(Tr) | (Tr<1 & T.o==0),
   T.h = []; T.v='';
end;

%end .. @tsp/private/tclear
