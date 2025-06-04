function r = poles(A,varargin);
%POLES   Poles of a right-den fraction
%
% The command
%   POLES(F)
% computes the poles of right-denominator fraction F,
% i.e. the roots of its denominator.
%
% The option  'all'  requires to compute the infinite poles
% as well. The options METHOD and TOL may specify the method
% of roots computation and the tolerance, as in POL/ROOTS.
%
% For fractions in 'z' or in 'z^-1', it is possible to prescribe
% the variable for computing poles, the default being  F.v .
% So, e.g.
%    POLES(z/(z-0.5))  or  POLES(z/(z-0.5),'z') or
%    POLES(1/(1-0.5*zi),'z')
% yields 0.5, whereas
%    POLES(1/(1-0.5*zi))  or  POLES(1/(1-0.5*zi),'zi') or
%    POLES(z/(z-0.5),'zi')  
% yields 2.
%
% See also POL/ROOTS.

%       Author: J.Jezek, 22-Nov-2002
%       Copyright(c) 2002 by Polyx, Ltd.
%       % Revision $  $ Date 28-Feb-2003  warning  $

global PGLOBAL;

if ~isa(A,'rdf'),
   error('Some argument but not 1st is invalidly rdf.');
end;

var = A.frac.v;
allr = '';
met = 'det';
tol = PGLOBAL.ZEROING;
recip = 0;

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
                           recip = 1; var = arg;
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

eval('A = coprime(A,tol);','error(peel(lasterr));');
if recip, A = reverse(A); A = coprime(A,tol);
end;

eval('r = roots(A.frac.den,var,allr,met,tol);','error(peel(lasterr));');

if ~isempty(allr),
   A.frac.v = 'z'; 
   A = reverse(A);
   A = reduce(A,tol);
   [ND,NL] = tdeg(A.frac.num,'col');
   [DD,DL] = tdeg(A.frac.den,'col');
   num = NL*diag(zi.^ND);
   den = DL*diag(zi.^DD);
   Azero = coprime(num/den);
   rzero = roots(Azero.frac.den);
   saved_w = warning; warning off;
   rinf = 1./rzero;
   warning(saved_w);
   r = [r;rinf];
end;

%end .. @rdf/poles         
