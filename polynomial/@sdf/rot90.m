function B = rot90(A,k)
%ROT90     Rotate the scalar-den fraction by 90 degrees
%           B = ROT90(A),    B = ROT90(A,K)
%
% ROT90(A) is the 90 degree counter-clockwise rotation of matrix A.
% ROT90(A,K) is the K*90 degree rotation of A. K is an integer.
%
% See also SDF/FLIPLR, SDF/FLIPUD.

%       Author:  J. Jezek  08-Feb-2000
%       Copyright (c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 26-Apr-2000 $
%                     $ Date 06-Nov-2000 $
%                     $ Date 30-Sep-2002 $
%                     $ Date 14-Oct-2002 $

global PGLOBAL;

if nargin==1, k = 1;
else
   if ~isa(k,'double'),
      error('Invalid 2nd argument.');
   end;
end;

B = 0;
eval ('B = sdf(rot90(A.frac.num,k),A.frac.den);','error(peel(lasterr));');

props(B,A.frac.c,A.frac.tc,A.frac.r,A.frac.p,A.frac.tp);
if strcmp(PGLOBAL.COPRIME,'cop'),
   B = coprime(B);
end;
if strcmp(PGLOBAL.REDUCE,'red'),
   B = reduce(B);
else
   B = smreduce(B);
end;
if strcmp(PGLOBAL.DEFRACT,'defr'),
   B = defract(B);
end;

%end .. @sdf/rot90
