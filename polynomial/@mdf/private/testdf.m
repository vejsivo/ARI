function [td,F,G] = testdf(F,G)
%TESTDF    Test dimensions of matrix-den fractions
%             [TD,F,G] = TESTDF(F,G)
% The command tests whether the dimensions of F,G are the same
% (resulting TD = 1). Admissible is also a case when one of the arguments
% is scalar. Its dimensions are changed according to the other argument,
% all entries being filled with the scalar.

%     Author: J. Jezek  03-Feb-2000
%     Copyright(c) 2000 by Polyx,Ltd.
%     $ Revision $  $ Date: 31-Jul-2000 $

global PGLOBAL;

savdefr = PGLOBAL.DEFRACT; PGLOBAL.DEFRACT = 'ndefr';
td = 1;
if ~all(F.frac.s==G.frac.s),
   if all(F.frac.s==1),
      F = repmat(F,G.frac.s);
   elseif all(G.frac.s==1),
      G = repmat(G,F.frac.s);
   else
      td = 0;
   end;
end;
PGLOBAL.DEFRACT = savdefr;

%end .. @mdf/private/testdf
