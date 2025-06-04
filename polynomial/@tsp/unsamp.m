function B = unsamp(A,T,varargin)
%UNSAMP    Unsample discrete time two-sided polynomial
%
% For two-sided polynomial A, the command
%      B = UNSAMP(A,T,TAU)  
% where scalars T,TAU are sampling period, T>=0, and
% sampling phase, 0<=TAU<T, always results in error.
% The argument TAU is optional, the default is zero.
%
% This macro exists only for completeness.
% See also FRAC/UNSAMPLE.

%        Author: J.Jezek  30-Nov-2000
%        Copyright(c) 2000 by Polyx, Ltd.

ni = nargin;
if ni<2,
   error('Not enough input arguments.');
end;
if ~isa(T,'double') | length(T)~=1 | ~isreal(T) | T<=0,
   error('Invalid sampling period.');
end;

tau = 0; lv = length(varargin);
for i = 1:lv,
   arg = varargin{i};
   if isa(arg,'double'),
      if ~isempty(arg),
         if ndims(arg)==2 & any(size(arg)==1) & ...
               isreal(arg) & all(arg>=0) & all(arg<T),
            tau = arg;
         else
            error(['Invalid ',nth(i+2),' argument.']);
         end;
      end;
   elseif isa(arg,'tsp'),
      error(['Invalid ',nth(i+2),' argument.']);
   end;
end;   

lv = length(varargin);
eval('B = unsamp(sdf(A),T,varargin{1:lv});', ...
   'error(peel(lasterr));');

%end .. tsp/unsamp

