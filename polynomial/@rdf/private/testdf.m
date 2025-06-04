function [td,F,G] = testdf(F,G)
%TESTDF    Test dimensions of right-den fractions
%             [TD,F,G] = TESTDF(F,G)
% The command tests whether the dimensions of fractions F,G are the same
% (resulting TD = 1). Admissible is also a case when one of the arguments
% is scalar. Its dimensions are changed according to the other argument,
% all entries being filled with the scalar. In other cases, TD = 0 .
%

%     Author: J. Jezek  12-May-2000
%     Copyright(c) 2000 by Polyx,Ltd.

global PGLOBAL;

savdefr = PGLOBAL.DEFRACT; PGLOBAL.DEFRACT = 'ndefr';

td = 1;
if ~all(F.frac.s==G.frac.s),
   if all(F.frac.s==1),
      F = repmat(F,G.frac.s);
   elseif all(Q.s==1),
      G = repmat(G,F.frac.s);
   else
      td = 0;
   end;
end;
PGLOBAL.DEFRACT = savdefr;

%end .. @rdf/private/testdf