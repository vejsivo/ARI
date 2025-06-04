function [b,blog] = integral(a,var)
%INTEGRAL     Integral of constant
%
% The command  B = INTEGRAL(A,VAR) for constant A
% (i.e. standard Matlab matrix), and for variable VAR,
% returns polynomial A*VAR. The second argument is
% optional. If present, it must be one of 's', 'p',
% 'z','z^-1','q','d'. The default is the standard
% global variable.
%
% This macro exists only for completeness.
% See also POL/INTEGRAL, TSP/INTEGRAL.

%         Author: J.Jezek, 10-Oct-2000
%         Copyright(c) 2000 by Polyx, Ltd.
%         $ Revision $  $ Date 22-Sep-2002 $

global PGLOBAL;
eval('PGLOBAL.VARIABLE;', 'painit;');

ni = nargin;
if ni<1,
   error('Not enough input arguments.');
end;
if ni==2,
   if isa(var,'char'),
      if isempty(var),
         var = PGLOBAL.VARIABLE;
      elseif isempty(strmatch(var,{'z';'z^-1';'zi';'s';'p';'q';'d'},'exact')),
         error('Invalid integration variable.');
      end;
      w = pol([0 1],1,var);
   else
      eval('w = pol(var);', ...
         'error(''Invalid integration variable.'');');
      [vs1,vs2,vd] = size(w);
      if all([vs1,vs2,vd]==1) & all(w.c(:,:)==[0,1]),
      else
         error('Invalid integration variable.');
      end;
   end;
else
   var = PGLOBAL.VARIABLE;
   w = pol([0 1],1,var);
end;

if ~((isa(a,'double') & ndims(a)<=2) | isa(a,'frac')),
   error('Invalid 1st argument.');
end;
notdouble = logical(0);
eval('a = double(a);', 'notdouble = logical(1);');
if notdouble,
   notpol = logical(0);
   eval('a = pol(a);', 'notpol = logical(1);');
   if notpol,
      nottsp = logical(0);
      eval('a = tsp(a);', 'nottsp = logical(1);');
      if nottsp,
         error('Invalid 1st argument.');
      end;
      if ~(strcmp(w.v,'z^-1')|strcmp(w.v,'z')),
         w.v = 'z';
         warning('Integration variable changed to ''z''.');
      end;
   end;
   if nargout==2,
      eval('[b,blog] = integral(a,w);', ...
         'error(peel(lasterr));');
   else
      eval('b = integral(a,w);', ...
         'error(peel(lasterr));');
   end;
else
   b = a*w;
   blog = pol(a*0);
end;

%end .. integral
