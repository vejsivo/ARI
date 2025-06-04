function B = chsamp(A,varargin)
%CHSAMP    Change sampling of discrete time polynomial
%
% For polynomial A, the command  B = CHSAMP(A,T,TN,TAU,TAUN)
% where scalars  T,TN,TAU,TAUN  are old and new sampling
% period and old and new sampling phase, in all nontrivial
% cases results in error.
%
% The arguments T,TN,TAU,TAUN are optional, they may be
% omitted or given by [].
%
% This macro exists only for completeness.
% See also RDF/CHSAMP, LDF/CHSAMP, MDF/CHSAMP, SDF/CHSAMP.

%        Author: J.Jezek  05-Feb-2001
%        Copyright(c) 2001 by Polyx, Ltd.

ni = nargin;
if ni<1,
   error('Not enough input arguments.');
end;
if ~isa(A,'pol'),
   error('Some argument but not 1st is invalidly polynomial.');
end;

lv = length(varargin);
eval('B = chsamp(sdf(A),varargin{1:lv});', ...
   'error(peel(lasterr));');
B = defract(B);

%end .. @pol/chsamp

