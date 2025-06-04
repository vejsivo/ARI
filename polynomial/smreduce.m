function G = smreduce(F)
%SMREDUCE    Small reduce constant
%
% The command  G = SMREDUCE(F) returns  G = F .
%
% This macro exists only for completeness.
% See also LDF/SMREDUCE, RDF/SMREDUCE, MDF/SMREDUCE, SDF/SMREDUCE.

%      Author:  J. Jezek  30-Sep-2002
%      Copyright(c) 2002 by Polyx, Ltd.

if nargin==0,
   error('Not enough input arguments.');
end;
if ~isa(F,'double') | ndims(F)>2,
   error('Invalid 1st argument.');
end;

G = F;

%end .. smreduce
