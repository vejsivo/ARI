function r = roots(A,varargin);
%ROOTS   Roots of a scalar-den fraction
%
% The command
%   ROOTS(F)
% computes the roots (zero-points) of scalar-denominator
% fraction F.
%
% The option  'all'  requires to compute the infinite roots
% as well. The options METHOD and TOL may specify the method
% of roots computation and the tolerance, as in POL/ROOTS.
%
% For fractions in 'z' or in 'z^-1', it is possible to prescribe
% the variable for computing roots, the default being F.v .
% So, e.g.
%    ROOTS((z-0.2)*inv(z-0.5))  or
%    ROOTS((z=0.2)*inv(z-0.5),'z') or
%    POLES((1-0.2*zi)*inv(1-0.5*zi),'z')
% yields 0.2, whereas
%    ROOTS((1=0.2*zi)*inv(1-0.5*zi))  or
%    ROOTS((1-0.2*zi)*inv(1-0.5*zi),'zi') or
%    ROOTS((z-0.2)*inv(z-0.5),'zi')  
% yields 5.
%
% See also POL/ROOTS.

%       Author: J.Jezek, 22-Nov-2002
%       Copyright(c) 2002 by Polyx, Ltd.

global PGLOBAL;

if ~isa(A,'sdf'),
   error('Some argument but not 1st is invalidly sdf.');
end;

var = A.frac.v;
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
               if strcmp(arg,'zi'), arg = 'z^-1';
               end;
               I = strmatch(arg,{'s';'p';'z^-1';'d';'z';'q'},'exact');
               if ~isempty(I),
                  if ~isempty(var),
                     if ~strcmp(var,arg),
                        if (strcmp(var,'z') & strcmp(arg,'z^-1')) | ...
                              (strcmp(var,'z^-1') & strcmp(arg,'z')),
                           var = arg;
                        else
                           error('Invalid variable for roots.');
                        end;
                     end;
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
if A.frac.s(1) >= A.frac.s(2),
   eval('A = ldf(A,tol);','error(peel(lasterr));');
else
   eval('A = rdf(A,tol);','error(peel(lasterr));');
end;

eval('r = roots(A,var,allr,met,tol);','error(peel(lasterr));');

%end .. @sdf/roots
