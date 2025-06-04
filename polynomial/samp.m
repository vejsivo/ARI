function B = samp(A,T,varargin)
%SAMP    Sample continuous time constant
%
% For constant A (i.e. the standard Matlab Matrix),
% the command   B = SAMP(A,T,TAU)  
% where scalars T,TAU are sampling period and
% sampling phase, T,TAU >= 0, in nontrivial cases
% results in error when TAU is zero, and returns 0
% when TAU is nonzero.
%
% The arguments T,TAU are optional, the defaults are
%  T = 1, TAU = 0 .
%
% This macro exists only for completeness.
% See also POL/SAMP, RDF/SAMP, LDF/SAMP,
% MDF/SAMP, SDF/SAMP.

%         Author: J. Jezek  19-Jun-2000
%         Copyright(c) by Polyx, Ltd.
%         $ Revision $  $ Date 17-Oct-2000 $
%                       $ Date 14-Dec-2000 $
%                       $ Date 26-Jan-2001 $
%                       $ Date 25-Jul-2002 $

ni = nargin;
if ni<1,
   error('Not enough input arguments.');
end;
if ~isa(A,'double') | ndims(A)~=2,
   error('Invalid 1st argument.');
end;

if ni<2, T = 1;
elseif ~isa(T,'double') | length(T)~=1 | ~isreal(T) | T<0,
   error('Invalid sampling period.');
end;

tau = 0; lv = length(varargin);
for i = 1:lv,
   arg = varargin{i};
   if isa(arg,'double'),
      if ~isempty(arg),
         if ndims(arg)==2 & any(size(arg)==1) & ...
               isreal(arg) & all(arg>=0),
            tau = arg;
         else
            error(['Invalid ',nth(i+2),' argument.']);
         end;
      end;
   end;
end;   

if ~isempty(A) & any(any(A~=0)) & any(tau==0),
   error('In time domain, sampling of Dirac impulse.');
end;

BB = zeros(size(A));
ltau = length(tau);
if ltau==1,
   B = BB;
else
   B = cell(size(tau));
   for k = 1:ltau,
      B{k} = BB;
   end;
end;

%end .. samp
