function Y = mvalue(F,X)
%MVALUE     Matrix value of constant
%
% For constant F (i.e. standard Matlab matrix),
% the command  Y = MVALUE(F,X)  returns F or
% F*eye matrix of dimension of X. Argument X
% plays no role but it must be constant square
% matrix. It may be also omitted.
%
% This macro exists only for completeness.
% See also POL/MVALUE, TSP/MVALUE, RDF/MVALUE,
% LDF/MVALUE, MDF/MVALUE, SDF/MVALUE.

%        Author: J.Jezek, 05-Oct-2000
%        Copyright(c) 2000 by Polyx, Ltd.
%        $ Revision $ Date 23-Apr-2003 $

ni = nargin;
if ni<1,
   error('Not enough input arguments.');
end;

if ~isa(F,'double') | ndims(F)>2,
   error('Invalid 1st argument.');
end;

if ni==2,
   eval('Y = mvalue(pol(F),X);','error(peel(lasterr));');
else
   eval('Y = mvalue(pol(F));','error(peel(lasterr));');
end;

%end .. mvalue

      