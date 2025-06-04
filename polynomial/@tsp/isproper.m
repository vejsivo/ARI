function t = isproper(F,varargin)
%ISPROPER     Test if two-sided polynomial is proper
%
% For two-sided polynomial F, the expression
%  ISPROPER(F)  returns 1 if F is proper, otherwise ).
% The properness means behaviour in point  z = Infinity.
% 
% An optional argument 'strictly' requires testing whether
% the fraction is strictly proper. It may be shortened down to 's'.
%  
% An optional argument TOL may specify the zeroing tolerance
% to be used instead of the standard one.

%     Author:  J. Jezek, 25-Sep-2002
%     Copyright(c) 2002 by Polyx, Ltd.

if ~isa(F,'tsp'),
   error('Invalid 1st argument.');
end;

lv = length(varargin);
eval('t = isproper(sdf(F),varargin{1:lv});', ...
   'error(peel(lasterr));');

%end .. @tsp/isproper
