function y = isproper(N,D,varargin)
%ISPROPER   Test if constant matrix is proper
%
% For constant matrices N,D, the commands
%    ISPROPER(N,D)
%    ISPROPER(N,D,'l') 
% return 1 if the left matrix fraction D\N is a proper rational 
% function, and 0 if it is not. The command
%    ISPROPER(N,D,'r') 
% returns 1 if the right matrix fraction N/D is a proper rational 
% function, and 0 otherwise. The command
%    ISPROPER(N,D,'strictly'), 
% possibly combined with the 'l' or 'r' option, tests if the matrix 
% fraction is strictly proper. The string 'strictly' may be shortened 
% to 'strict', 'str' or 's'.
%
% If D is omitted or empty then D = 1  is taken.
%
% A constant matrix is always taken as proper.
% It is taken as strictly proper only when it is empty or zero.
%
% This macro exists only for completeness.
% See also POL/ISPROPER, RDF/ISPROPER, RDF/ISPROPER, FRAC/ISPROPER.

%      Author: J. Jezek, 10-Jul-2001
%      Copyright(c) 2001 by Polyx, Ltd.
%      $ Revision $  $ Date 25-Jul-2002 $

opt = 'l';
strict = 0;
twoarg = 0;

ni = nargin;
if ni<1,
   error('Not enough input arguments.');
end;
eval('N = double(N);', 'error(peel(lasterr));');

if ni>=2 & ~isempty(D),
   if ischar(D),
      if D(1)=='s',
         strict = 1;
      else
         opt = D;
      end;
   else
      eval('D = double(D);', 'error(peel(lasterr));');
      twoarg = 1;
   end;
end;

lv = length(varargin);
for i = 1:lv,
   arg = varargin{i};
   if isnumeric(arg) & all(size(arg))==1 & ...
         isreal(arg) & arg>=0 & arg<=1,
      tol = arg;
   elseif ischar(arg),
      if arg(1)=='s',
         strict = 1;
      else
         opt = arg;
      end;
   else
      error(['Invalid ',nth(i+2),' argument.']);
   end;
end;

if ~twoarg, D = 1;
end;

% size consistency:
if ~(strcmp(opt,'l') | strcmp(opt,'r')),
   error('Invalid command option.');
end;
eval('[N,D] = testdnd(N,D,opt);', ...
   'error(peel(lasterr));');

% D must have full rank:
if rank(D) < size(D,1),
   warning('Denominator is singular.');
   y = logical(0);
   return;
end;
    
if strict,
   if isempty(N) | all(N==0),
      y = logical(1);
   else
      y = logical(0);
   end;
else
   y = logical(1);
end;

%end .. isproper
