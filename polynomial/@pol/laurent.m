function [Q,R] = laurent(F,varargin);
%LAURENT      Laurent series of polynomial
%
% For polynomial F in variables 'z', 'q', 's' or 'p',
% the command  Q = LAURENT(F,DEGREE)  returns  Q = F,
% DEGREE being ignored.
%
% For polynomial in 'z^-1' or 'd',
% the command  Q = LAURENT(F,DEGREE)  returns  Q = F,
% possibly truncated to DEGREE.
%
% For polynomial in 's' or 'p', the command
%   [Q,R] = LAURENT(F,DEGREE),  returns  Q = positive 
% powers of F, and  R ... three-dimensional array,
%  R(:,:,1) = F{0},  R(:,:,I) = 0  for I = 2,...DEGREE+1.
%
% The default for DEGREE is  DEG(F).
%
% This macro exists only for completeness.
% See also RDF/LAURENT, LDF/LAURENT.

%      Author: J.Jezek, 18-Jul-2000
%      Copyright(c) 2000 by Polyx, Ltd.
%      $ Revision $  $ Date 24-Sep-2002  comments $

if nargin<1,
   error('Not enough input arguments.');
end;
eval('F = pol(F);','error(peel(lasterr));');

lv = length(varargin);
if nargout<=1,
   eval('Q = laurent(sdf(F),varargin{1:lv});', ...
      'error(peel(lasterr));');
else
   eval('[Q,R] = laurent(sdf(F),varargin{1:lv});', ...
      'error(peel(lasterr));');
end;

%end .. @pol/laurent


   