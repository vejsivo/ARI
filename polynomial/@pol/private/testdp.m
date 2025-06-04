function [td,P,Q] = testdp(P,Q)
%TESTDP    Test dimensions of polynomials
%             [TD,P,Q] = TESTDP(P,Q)
% The command tests whether the dimensions of polynomials P,Q are the same
% (resulting TD = 1). Admissible is also a case when one of the arguments
% is scalar. Its dimensions are changed according to the other argument,
% all entries being filled with the scalar. In other cases, TD = 0 .

%     Author: J. Jezek  12-May-2000
%     Copyright(c) 2000 by Polyx,Ltd.

td = 1;
if ~all(P.s==Q.s),
   if all(P.s==1),
      P = repmat(P,Q.s);
   elseif all(Q.s==1),
      Q = repmat(Q,P.s);
   else
      td = 0;
   end;
end;

%end .. @pol/private/testdp
