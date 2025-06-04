function r = poles(A,varargin);
%POLES   Poles of a left-den fraction
%
% The command
%   POLES(F)
% computes the poles of left-denominator fraction F,
% i.e. the roots of its denominator.
%
% The option  'all'  requires to compute the infinite poles
% as well. The options METHOD and TOL may specify the method
% of roots computation and the tolerance, as in POL/ROOTS.
%
% For fractions in 'z' or in 'z^-1', it is possible to prescribe
% the variable for computing poles, the default being F.v .
% So, e.g.
%    POLES((z-0.5)\z)  or  POLES((z-0.5)\z,'z') or
%    POLES((1-0.5*zi)\1,'z')
% yields 0.5, whereas
%    POLES((1-0.5*zi)\1)  or  POLES((1-0.5*zi)\1,'zi') or
%    POLES((z-0.5)\z,'zi')  
% yields 2.
%
% See also POL/ROOTS.

%       Author: J.Jezek, 22-Nov-2002
%       Copyrigh(c) 2002 by Polyx, Ltd.
%       $ Revision $  $ Date 28-Feb-2003  warnings  $

global PGLOBAL;

if ~isa(A,'ldf'),
   error('Some argument but not 1st is invalidly ldf.');
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
   [ND,NL] = tdeg(A.frac.num,'row');
   [DD,DL] = tdeg(A.frac.den,'row');
   num = diag(zi.^ND)*NL;
   den = diag(zi.^DD)*DL;
   Azero = coprime(den\num);
   rzero = roots(Azero.frac.den);
   saved_w = warning; warning off;
   rinf = 1./rzero;
   warning(saved_w);
   r = [r;rinf];
end;

%end .. @ldf/poles         
