function B = dilate(A,varargin)
%DILATE     Dilate discrete-time constant
%
% For constant A (i.e. the standard Matlab matrix), the command
%    B = DILATE(A,K,H)
% where scalar integers K,H are dilating period,  K>=1,  and
% dilating phase,  0<=H<K,  returns constant, polynomial or
% fraction  B(z^-1) = z^-H * A. The argument H is optional,
% the default value being 0.
%
% The variable symbol of B is the same as that of A. If A has
% no symbol then the symbol of B may be taken from an optional
% (fourth) argument VAR or from the global discrete variable
% symbol.
%
% This macro exists only for completeness.
% See also POL/DILATE, TSP/DILATE, FRAC/DILATE, RESAMP.

%        Author: J.Jezek, 05-Oct-2000
%        Copyright(c) 2000 by Polyx, Ltd.
%        $ Revision $  $ Date 28-Feb-2003 $

ni = nargin;
if ni<1,
   error('Not enough input arguments.');
end;
if ~isa(A,'double') | ndims(A)~=2,
   error('Invalid 1st argument.');
end;

eval('B = dilate(pol(A),varargin{:});', ...
        'error(peel(lasterr));');

%end .. dilate
