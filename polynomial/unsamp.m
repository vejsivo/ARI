function B = unsamp(A,varargin)
%UNSAMP    Unsample discrete time constant
%
% For constant A (i.e. the standard Matlab matrix),
% the command   B = UNSAMP(A,T,TAU)  
% where scalars T,TAU are sampling period and
% sampling phase, T,TAU >= 0, in all nontrivial cases
% results in error.
%
% The arguments T, TAU are optional, they may be
% omitted or given by [].
%
% This macro exists only for completeness.
% See also FRAC/UNSAMP.

%        Author: J.Jezek  17-Oct-2000
%        Copyright(c) 2000 by Polyx, Ltd.
%        $ Revision $  $ Date 14-Dec-2000 $
%                      $ Date 27-Dec-2000 $

ni = nargin;
if ni<1,
   error('Not enough input arguments.');
end;
if ~isa(A,'double') | ndims(A)~=2,
   error('Invalid first argument.');
end;

lv = length(varargin);
if lv>=1,
   arg = varargin{1};
   if ~isempty(arg),
      if ~isa(arg,'double') | length(arg)~=1 | ~isreal(arg) | arg<0,
         error('Invalid sampling period.');
      end;
   end;

   for i = 2:lv,
      arg = varargin{i};
      if isa(arg,'double'),
         if ~isempty(arg),
            if length(arg)==1 & isreal(arg) & arg>=0,
            else
               error(['Invalid ',nth(i+2),' argument.']);
            end;
         end;
      end;
   end;
end;


if isempty(A) | all(all(A==0)),
   B = A; return;
end;

error('Invalid first argument; has not required zero point.');

%end .. unsamp

