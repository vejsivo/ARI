function Z = null(T,varargin)
%NULL  Nullspace of two-sided polynomial
%
% The command
%    Z = NULL(F)
% computes a polynomial basis for the right null space of the
% two-sided polynomial T. The result Z is a polynomial matrix.
% If F has full column rank then Z is an empty polynomial matrix.
%
% Optional input arguments DEGREE or TOL may be given as in POL/NULL.
%
% See also POL/NULL.

%      Author:  J. Jezek, 13-Apr-2000
%      Copyright(c) 2000 by Polyx, Ltd.

if ~isa(T,'tsp'),
   error('Invalid arguments.');
end;

ni = length(varargin);
eval('Z = null(T.p,varargin{1:ni});', ...
   'error(peel(lasterr));');

%end .. @tsp/null
