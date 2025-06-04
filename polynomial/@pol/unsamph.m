function  F = unsamph(G,T,varargin)
%UNSAMPH     Un- sample and hold discrete time polynomial
%
% Polynomial G is supposed to be a result of sampling
% and holding of continuous time polynomial F with
% period T and phase TAU. The command
%     F = UNSAMPH(G,T,TAU)
% restores the original polynomial F.
%
% The only admissible possibility for G is to be constant.
% In this case, F = G, for all T and TAU. All other cases
% result in error. Argument T is compulsory, TAU is optional.
%
% This macro exists only for completeness.
% See also  FRAC/UNSAMPH.

%       Author: J.Jezek  01-Dec-2000
%       Copyroght(c) 2000 by Polyx, Ltd.

if nargin<2,
   error('Not enough input arguments.');
end;
if ~isa(G,'pol'),
   error('Invalid 1st argument.');
end;

lv = length(varargin);
eval('F = unsamph(sdf(G),T,varargin{1:lv});', ...
   'error(peel(lasterr));');
F = defract(F);

%end .. @pol/unsamph
