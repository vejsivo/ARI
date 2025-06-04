function r = poles(A,varargin);
%POLES   Poles of a polynomial matrix
%
% The command
%   POLES(P)
% computes the poles of polynomial matrix P. The result
% is empty, as only finite poles are computed. However,
% the option  'all'  requires to compute the infinite
% poles as well. 
%
% For polynomials in 'z' or 'z^-1', it is possible to prescribe
% the variable for computing poles, the default being P.v .
% So, e.g.
%   POLES(1+z)            yields  empty,
%   POLES(1+z,'all')      yields  Inf,
%   POLES(1+z,'zi')       yields  0,
%   POLES(1+z,'zi','all') yields  0.
%
% See also POL/ROOTS, RDF/POLES.

%       Author: J.Jezek, 22-Nov-2002
%       Copyright(c) 2002 by Polyx, Ltd.

global PGLOBAL;

if ~isa(A,'pol'),
   error('Some argument but not 1st is invalidly pol.');
end;

var = A.v;
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

if A.s(1) >= A.s(2),
   A = ldf(A);
else
   A = rdf(A);
end;

eval('r = poles(A,var,allr,met,tol);','error(peel(lasterr));');

%end .. @pol/poles         
