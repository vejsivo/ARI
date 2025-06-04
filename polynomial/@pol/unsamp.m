function B = unsamp(A,varargin)
%UNSAMP    Unsample discrete time polynomial
%
% For polynomial A, the command   B = UNSAMP(A,T,TAU)  
% where scalars T,TAU are sampling period and
% sampling phase, T,TAU >= 0, in all nontrivial cases
% results in error.
%
% The arguments T,TAU are optional, they may be omitted
% or given by [].
%
% This macro exists only for completeness.
% See also FRAC/UNSAMP.

%        Author: J.Jezek  30-Nov-2000
%        Copyright(c) 2000 by Polyx, Ltd.
%        $ Revision $  $ Date 14-Dec-2000 $
%                      $ Date 27-Dec-2000 $

ni = nargin;
if ni<1,
   error('Not enough input arguments.');
end;
if ~isa(A,'pol'),
   error('Some argument but not 1st is invalidly polynomial.');
end;

lv = length(varargin);
eval('B = unsamp(sdf(A),varargin{1:lv});', ...
   'error(peel(lasterr));');
B = defract(B);

%end .. @pol/unsamp

