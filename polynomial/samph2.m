function G = samph2(F,varargin)
%SAMPH2     Sample and second order hold continuous time constant
%
% For constant F, the command
%    G = SAMPH2(F,T,TAU,L)   returns constant, polynomial
% or fraction G, the result of second order holding with
% interval L and sampling with period T and phase TAU.
%
% For more details, see POL/SAMPH2, FRAC/SAMPH2.
%
%    Author: J.Jezek, 13-Jul-2000
%    Copyright(c) 2000 by Polyx,Ltd.
%    $ Revision $  $ Date 04-Oct-2000 $
%                  $ Date 01-Dec-2000 $
%                  $ Date 25-Jul-2002 $

if nargin<1,
   error('Not enough input arguments.');
end;
if ~isa(F,'double') | ndims(F)~=2,
   error('Invalid 1st argument.');
end;

lv = length(varargin);
eval('G = samph2(pol(F),varargin{1:lv});',...
   'error(peel(lasterr));');
if ~isa(G,'cell'),
   G = defract(G);
end;

%end .. samph2
