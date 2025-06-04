function B = resamp(A,k,h)
%RESAMP     Resample discrete time constant
%
% For constant A (i.e. the standard Matlab matrix), the command
%   B = RESAMP(A,K,H)
% where scalar integers K,H are resampling period, K>=1, and
% resampling phase, H>=0, (defaults being  K=2, H=0 ), returns A
% when H is zero, and returns 0 when H is nonzero. 
%
% H can also be, instead of scalar, a vector. In such
% a case, the result is a cell vector.
%
% This macro exists only for completeness.
%
% See also POL/RESAMP, TSP/RESAMP, RDF/RESAMP, LDF/RESAMP,
% MDF/RESAMP, RDF/RESAMP, DILATE.

%       Author:  J. Jezek 19-Jun-2000
%       Copyright(c) by 2000 Polyx, Ltd.
%       $ Revision $  $ Date 04-Oct-2000 $
%                     $ Date 01-Feb-2001 $
%                     $ Date 25-Jul-2002 $
%                     $ Date 28-Feb-2003 $

ni = nargin;
if ni<1,
   error('Not enough input arguments.');
end;
if ~isa(A,'double') | ndims(A)~=2,
   error('Invalid 1st argument.');
end;

if ni<2 | isempty(k),
   k = 2;
else
   if ~isa(k,'double') | length(k)~=1 | ~isreal(k) | ...
         floor(k)~=k | k<=0,
      error('Invalid resampling period.');
   end;
end;

if ni<3 | isempty(h),
   h = 0;
else
   if ~isa(h,'double') | ndims(h)~=2 | ~any(size(h)==1) | ...
         ~isreal(h) | any(floor(h)~=h) | any(h<0),
      error('Invalid resampling phase.');
   end;
end;

lh = length(h); sA = size(A);
if lh==1,
   if h==0, B = A;
   else B = zeros(sA);
   end;
else
   B = cell(size(h));
   for i = 1:lh,
      if h(i)==0, B{i} = A;
      else B{i} = zeros(sA);
      end;
   end;
end;

%end .. resamp
