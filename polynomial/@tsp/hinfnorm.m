function n = hinfnorm(A,varargin)
%HINFNORM    H-infinity norm of a two-sided polynomial
% 
% The command  HINFNORM(A)  computes H-infinity norm
% of two-sided polynomial A.
%
% See also POL/HINFNORM.

%      Author: J.Jezek 02-Jan-2002
%      Copyright(c) 2002 by Polyx, Ltd.

if ~isa(A,'tsp'),
    error('Some argument but not 1st is invalidly tsp.');
end;
eval('n = hinfnorm(pol(shift(A,tdeg(A))),varargin{:});',...
    'error(peel(lasterr));');

%end .. @tsp/hinfnorm
