function n = h2norm(A,varargin)
%H2NORM    H2-norm of a two-sided polynomial
% 
% The command  H2NORM(A)  computes H2-norm
% of two-sided polynomial A.
%
% See also POL/H2NORM.

%      Author: J.Jezek 02-Jan-2002
%      Copyright(c) 2002 by Polyx, Ltd.

if ~isa(A,'tsp'),
    error('Some argument but not 1st is invalidly tsp.');
end;
eval('n = h2norm(pol(shift(A,tdeg(A))),varargin{:});',...
    'error(peel(lasterr));');

%end .. @tsp/h2norm
