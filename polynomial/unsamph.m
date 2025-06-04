function  F = unsamph(G,T,varargin)
%UNSAMPH     Un- sample and hold discrete time constant
%
% Constant G (standard MATLAB matrix) is supposed to be
% a result of sampling and holding of continuous time
% constant F with period T and phase TAU. The command
%     F = UNSAMPH(G,T,TAU)
% restores the original polynomial F.
%
% It is  F = G,  for all T and TAU.
% Argument T is compulsory, TAU is optional.
%
% This macro exists only for completeness.
% See also  FRAC/UNSAMPH.

%       Author: J.Jezek  01-Dec-2000
%       Copyright(c) 2000 by Polyx, Ltd.

if nargin<2,
   error('Not enough input arguments.');
end;
if ~isa(G,'double'),
   error('Invalid first argument.');
end;

lv = length(varargin);
eval('F = unsamph(sdf(G),T,varargin{1:lv});', ...
   'error(peel(lasterr));');
F = declass(F);

%end .. unsamph
