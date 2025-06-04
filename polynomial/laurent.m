function [Q,R] = laurent(F,degree,varargin);
%LAURENT      Laurent series of a constant
%
% For a constant matrix F,
% the command  Q = LAURENT(F,DEGREE)  returns  Q = F,
% DEGREE being ignored.
%
% The command  [Q,R] = LAURENT(F,DEGREE)  returns  Q = 0,
% and  R ... three-dimensional array,  R(:,:,1) = F,
% R(:,:,I) = 0  for  I = 1,2,...DEGREE+1.
%
% The default for DEGREE is 0.
%
% This macro exist only for completeness.
% See also RDF/LAURENT, LDF/LAURENT, POL/LAURENT.

%       Author: J.Jezek, 18-Jul-2000
%       Copyright(c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 24-Sep-2002 $

if nargin<1,
   error('Not enough input arguments.');
end;
if ~isa(F,'double') | ndims(F)>2,
   error('Invalid 1st argument.');
end;

if nargin<2, degree = 0;
elseif ~isa(degree,'double') | length(degree)~=1 | ~isfinite(degree) |...
      ~isreal(degree) | floor(degree)~=degree | degree<0,
   error('Invalid degree.');
end;

if nargout<=1, Q = F;
else
   [sF1,sF2] = size(F);
   Q = zeros(sF1,sF2);
   R = zeros(sF1,sF2,degree+1);
   R(:,:,1) = F;
end;

%end .. laurent

   
   
   
   
   
