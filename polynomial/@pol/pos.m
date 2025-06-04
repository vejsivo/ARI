function PP = pos(P)
%POS  Positive part of a polynomial
%
% The command
%     P = POS(PP)
% for the polynomial PP, returns the part with positive
% powers only. It is always zero if PP is polynomial in z^-1.
%
% See also POL/NEG, POL/NPOS, POL/NNEG.

%        Author:  J.Jezek  11-8-99
%        Copyright (c) 1999 by Polyx, Ltd.

PP = P;
if ~isempty(P.d),
   Pv = P.v;
   if isempty(Pv) | strcmp(Pv,'z^-1'),
      Ps = P.s; PP = pol(zeros(Ps(1),Ps(2)));
   else
      U = P.c; U(:,:,1) = 0; PP.c = U;
   end;
end;

%end .. @pol/pos
