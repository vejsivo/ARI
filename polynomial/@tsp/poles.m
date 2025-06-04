function r = poles(A,varargin);
%POLES   Poles of a two-sided polynomial matrix
%
% The command
%   POLES(T)
% computes the poles of two-sided polynomial matrix T. The
% result is empty or contains 0, as only finite poles are
% computed. However, the option  'all'  requires to compute
% the infinite poles as well. 
%
% It is possible to prescribe the variable for computing
% poles, the default being 'z'. So, e.g.
%   POLES(z+z^-1)            yields  0,
%   POLES(z+z^-1,'all')      yields  0,Inf,
%   POLES(z+z^-1,'zi')       yields  0,
%   POLES(z+z^-1,'zi','all') yields  0,Inf.
%
% See also TSP/ROOTS, RDF/POLES.

%       Author: J.Jezek, 22-Nov-2002
%       Copyright(c) 2002 by Polyx, Ltd.

global PGLOBAL;

if ~isa(A,'tsp'),
   error('Some argument but not 1st is invalidly tsp.');
end;

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
               if strcmp(arg,'zi'), arg = 'z^-1';
               end;
               I = strmatch(arg,{'s';'p';'z^-1';'d';'z';'q'},'exact');
               if ~isempty(I),
                  if ~strcmp(var,arg),
                     if strcmp(arg,'z^-1'),
                        var = arg;
                     else
                        error('Invalid variable for roots.');
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

if A.s(1) >= A.s(2),
   A = ldf(A);
else
   A = rdf(A);
end;

eval('r = poles(A,var,allr,met,tol);','error(peel(lasterr));');

%end .. @tsp/poles         
