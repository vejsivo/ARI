function [x,y] = axyab(a,bl,bh,varargin)
%AXYAB   Non-symmetric equation solver
%
% The command
%    [X,Y] = AXYAB(A,BL,BH)
% solves the scalar equation
%    A'X + Y'A = BL' + BH
% with scalars A,BL,BH. The resulting X,Y are also
% scalars. The solution is not unique.
%
% This macro exist only for completeness. For more
% details and for possible further arguments,
% see POL/AXYAB.

%     Author: J.Jezek 06-Aug-2002
%     Copyright(c) 2001 by Polyx, Ltd.

ni = nargin;
if ni<3,
   error('Not enough input arguments.');
elseif ni==3,
   eval('[x,y] = axyab(pol(a),pol(bl),pol(bh),0);', ...
      'error(peel(lasterr));');
else
   eval('[x,y] = axyab(pol(a),pol(bl),pol(bh),varargin{:});', ...
      'error(peel(lasterr));');
end;

%end .. axyab
