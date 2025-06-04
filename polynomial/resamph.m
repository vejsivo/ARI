function G = resamph(F,varargin)
%RESAMPH     Resample and hold discrete time constant
%
% For constant F, the command
%  G = RESAMPH(F,K,H,L)  returns polynomial or fraction G,
% the result of holding with interval L and resampling with
% period K and phase H.
%
% For more details, see POL/RESAMPH, FRAC/RESAMPH.

%    Author: J.Jezek, 13-Jun-2000
%    Copyright(c) 2000 by Polyx,Ltd.
%    $ Revision $  $ Date 04-Oct-2000 $
%                  $ Date 25-Jul-2002 $

if nargin<1,
   error('Not enough input arguments.');
end;
if ~isa(F,'double') | ndims(F)~=2,
   error('Invalid 1st argument.');
end;

lv = length(varargin);
eval('G = resamph(pol(F),varargin{1:lv});',...
   'error(peel(lasterr));');

%end .. resamph
