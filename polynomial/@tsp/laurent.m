function [Q,R] = laurent(F,varargin);
%LAURENT      Laurent series of two-sided polynomial
%
% For two-sided polynomial F, the command
%   Q = LAURENT(F,DEGREE)  returns polynomial or
% two-sided polynomial  Q = F,  possibly truncated
% to DEGREE in negative powers of 'z'.
% 
% The default for DEGREE is  -TDEG(F).
%
% This macro exists only for completeness.
% See also RDF/LAURENT, LDF/LAURENT, POL/LAURENT.

%       Author: J.Jzek, 24-Sep-2002
%       Copyright(c) 2002 by Polyx, Ltd.

if nargin<1,
   error('Not enough input arguments.');
end;
eval('F = tsp(F);','error(peel(lasterr));');

lv = length(varargin);
if nargout<=1,
   eval('Q = laurent(sdf(F),varargin{1:lv});', ...
      'error(peel(lasterr));');
else
   eval('[Q,R] = laurent(sdf(F),varargin{1:lv});', ...
      'error(peel(lasterr));');
end;

%end .. @tsp/laurent
