function n = numel(varargin)
%NUMEL  Overloaded function for use with POL object.
%
% NUMEL(P) In the case of a general overloaded SUBSREF function, NUMEL is
%          used to compute the number of expected outputs (NARGOUT) returned
%          from SUBSREF. In the case of POL objects, there is always ONLY ONE
%          output argument, so NUMEL returns 1. For example, if P is a scalar
%          polynomial, then A = P{3:5} returns a single (!!!) vector composed
%          of coefficients with powers 3 to 5. Without overloading, POL/SUBSREF
%          would expect 3 output arguments, which would be conflicting.
%
%          For a general overloaded SUBSASGN function, NUMEL is used to compute
%          the number of expected inputs (NARGIN) to be assigned using SUBSASGN.
%          For POL objects, this overloaded NUMEL should also return 1. For example,
%          A{2:4} = [0 0 0] makes the coefficients with the powers of 2 to 4 zero.
%          Note, that the vector on the right side is again a single input argument,
%          so POL/SUBSASGN must also expect only one input. Without overloading,
%          it would expect three input parameters, which would be conflicting.

% Author(s): Z.Hurak, M.Hromcik        Oct.4 2001.
% Copyright (c) 2001 by PolyX, Ltd.
% Revision: 1.0.0 $Date: Oct.4 2001$

n = 1;

%end .. @pol/numel
