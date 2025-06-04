function L = llm(varargin)
%LLM   Least left multiple
%
% The command
%    L = LLM(N1, N2, .., Nk) 
% computes a least left common multiple L of several polynomial 
% matrices N1, N2, .., Nk (k > 0) that all have the same numbers
% of columns.
%
% Matrices Mi such that L = Mi*Ni may be recovered with Mi = XAB(Ni,L).
%
% A tolerance TOL may be specified as an additional input argument.
%
% See also: LRM.

%     Author: D. Henrion,  May 26, 1998.
%     Updated to 3.0 by D. Henrion, August 30, 2000.
%     Copyright 1998-2000 by Polyx, Ltd.
%     Modified by J.Jezek 24-Jul-2001, error check

% The function is dual to LRM.

% Transpose input arguments.

argin = {};
if nargin>0,
 for i = 1:nargin,
  argin{i} = varargin{i}.';
 end;
end;

% Call LRM.

eval('L = lrm(argin).'';', 'error(peel(lasterr));');

%end .. llm

