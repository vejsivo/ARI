function B = rot90(A,k)
%ROT90     Rotate right-den fraction by 90 degrees
%           B = rot90(A),    B = rot90(A,k)
%
% ROT90(A) is the 90 degree counter-clockwise rotation of matrix A.
% ROT90(A,K) is the K*90 degree rotation of A. K is an integer.
%
% See also RDF/FLIPLR, RDF/FLIPUD.

%       Author:  J. Jezek  08-Feb-2000
%       Copyright (c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 25-Apr-2000 $
%                     $ Date 01-Nov-2000 $
%                     $ Date 30-Sep-2002 $

global PGLOBAL;

if nargin==1, k = 1;
else
   if ~isa(k,'double'),
      error('Invalid 2nd argument.');
   end;
end;

A = sdf(A); B = 0;
eval('B = rot90(A,k);','error(peel(lasterr));');
B = rdf(B);

props(B,A.p,A.tp);
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

%end .. @rdf/rot90
