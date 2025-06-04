function [H,Q] = laplace(F,T,varargin)
%LAPLACE    Inverse Laplace transform of constant or polynomial
%
% For F constant or polynomial in variable 's' or 'p',
% the command    [H,Q] = LAPLACE(F,T,TAU,DEGREE)
% returns H=0, Q=F.
%
% This macro exists only for completeness.
% See also FRAC/LAPLACE.

%     Author: J.Jezek, 14-Jul-2000
%     Copyright(c) 2000 by Polyx, Ltd.

if nargin<2,
   error('Not enough input arguments.');
end;
lv = length(varargin);
eval('[H,Q] = laplace(sdf(F),T,varargin{1:lv});', ...
   'error(peel(lasterr));');

%end .. laplace
