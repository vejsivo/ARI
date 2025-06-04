function Y = value(F,X)
%VALUE     Value of constant
%
% For constant F (i.e. standard Matlab matrix),
% the command  Y = VALUE(F,X)  returns F or
% F * ones-matrix of dimension of X. Matrices
% F,X must have the same dimensions unless
% one of them is scalar. Matrix X may be omitted.
%
% This macro exists only for completeness.
% See also POL/VALUE, TSP/VALUE, RDF/VALUE,
% LDF/VALUE, MDF/VALUE, SDF/VALUE.

%        Author: J.Jezek, 05-Oct-2000
%        Copyright(c) 2000 by Polyx, Ltd.
%        $ Revision $  $ Date 27-Apr-2003 $

ni = nargin;
if ni<1,
   error('Not enough input arguments.');
end;
if ~isa(F,'double') | ndims(F)>2,
   error('Invalid 1st argument.');
end;
if ni==2,
   if ~isa(X,'double') | ndims(X)~=2,
      error('Invalid 2nd argument.');
   end;
   sF = size(F); sX = size(X);
   if all(sF==1),
      Y = repmat(F,sX);
   elseif any(sF~=sX) & any(sX~=1),
      error('Matrices not of the same dimensions.');
   else
      Y = F;
   end;
else
   Y = F;
end;   

%end .. value

      