function Z = null(F,varargin)
%NULL  Nullspace of right-den fraction
%
% The command
%    Z = NULL(F)
% computes a polynomial basis for the right null space of the
% right-den fraction F. The result Z is a polynomial matrix.
% If F has full column rank then Z is an empty polynomial matrix.
%
% Optional input arguments DEGREE or TOL may be given as in POL/NULL.
%
% See also POL/NULL.

%      Author:  J. Jezek, 13-Apr-2000
%      Copyright(c) 2000 by Polyx, Ltd.
%      $ Revision $  $ Date 25-Apr-2000 $
%                    $ Date 14-Oct-2002 $

if ~isa(F,'rdf'),
   error('Invalid arguments.');
end;

ni = length(varargin);
eval('Z = null(F.frac.num,varargin{1:ni});', ...
   'error(peel(lasterr));');

%end .. @rdf/null
