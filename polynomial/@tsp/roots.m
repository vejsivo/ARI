function r = roots(A,varargin);
%ROOTS  Roots of a two-sided polynomial matrix
%
% The command
%    ROOTS(T)
% computes the roots of a two-sided polynomial matrix T.
%
% The option  'all'  requires to compute the infinite roots
% as well. The options METHOD and TOL may specify the method
% of computation and the tolerance, as in POL/POOTS.
%
% The options 'z' or 'z^-1' may specify the variable for
% computing roots. The default is 'z'. So, e.g.
%    ROOTS(z - 4*z^-1)  or  ROOTS(z - 4*z^-1, 'z')
% yields  2,-2 ,  whereas
%    ROOTS(z - 4*z^-1, 'z^-1)  yields  0.5,-0.5 .
%
% See also POL/ROOTS.

%      Author: J.Jezek, 22-Nov-2002
%      Copyright(c) 2002 by Polyx, Ltd.

global PGLOBAL;

eval('A = tsp(A);','error(peel(lasterr));');
var = 'z';
allr = '';
met = 'det';
tol = PGLOBAL.ZEROING;

ni = nargin;
if ni>=2,
   for i = 2:ni,
      arg = varargin{i-1};
      if ~isempty(arg),
         if isa(arg,'pol'),
            [vs1,vs2,vd] = size(arg);
            if all([vs1,vs2,vd]==1) & all(arg.c(:,:)==[0 1]),
               arg = arg.v;
            else
               error(['Invalid ',nth(i),' argument.']);
            end;
         end;
         if ischar(arg),
            if strcmp(arg,'all'),
               allr = 'all';
            else
               I = strmatch(arg,{'s';'p';'z^-1';'d';'z';'q'},'exact');
               if ~isempty(I),
                  if strcmp(arg,'z'),
                  elseif strcmp(arg,'zi') | strcmp(arg,'z^-1'),
                     var = 'z^-1';
                  else
                     error('Invalid variable for roots.');
                  end;
               else
                  met = arg;
               end;
            end;
         elseif isnumeric(arg),
            tol = arg;
         else
            error(['Invalid ',nth(i),' argument.']);
         end;
      end;
   end;
end;

eval('r = roots(A.p,var,allr,met,tol);','error(peel(lasterr));');

%end .. @tsp/roots

