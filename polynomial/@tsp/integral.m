function [b,blog] = integral(a,var)
%INTEGRAL      Integral of two-sided polynomial
%
% The command  B = INTEGRAL(A,VAR) returns the integral
% of two-sided polynomial A. The integration constant
% is selected such that the absolute term of B is zero.
%
% The optional input argument VAR is the integration
% variable. If present, it must be one of 'z', 'z^-1'.
% The default is 'z'.
%
% The integral may contain a logarithmic term; its
% coefficient can be obtained in the second output
% argument.
%
% See also POL/INTEGRAL, TSP/DERIV.

%        Author: J.Jezek, 10-Oct-2000
%        Copyright(c) by Polyx, Ltd.
%        $ Revision $  $ Date 01-Aug-2001 sampl per $
%                      $ Date 30-Jun-2002 log term  $
%                      $ Date 22-Sep-2002           $

ni = nargin;
if ni<1,
   error('Not enough input arguments.');
end;
eval('a = tsp(a);', ...
   'error(''Invalid 1st argument.'');');
h = [];

if ni==2,
   if isa(var,'char'),
      if isempty(var),
         var = 'z';
      elseif isempty(strmatch(var,{'z';'z^-1';'zi'},'exact')),
         errot('Invalid integration variable.');
      end;
   else
      eval('var = pol(var);', ...
         'error(''Invalid integration variable.'');');
      [vs1,vs2,vd] = size(var);
      if all([vs1,vs2,vd]==1) & all(var.c(:,:)==[0,1]),
         h = var.h; var = var.v;
         if isempty(strmatch(var,{'z';'z^-1';'zi'},'exact')),
            error('Invalid integration variable.');
         end;
      else
         error('Invalid integration variable.');
      end;
   end;
else
   var = 'z';
end;

if ~isempty(h) & isfinite(h) & ~isempty(a.h) & isfinite(a.h),
   if a.h~=h,
      warning('Inconsistent sampling periods.');
      a.h = NaN;
   end;
end;
   
if isempty(a.h) & ~isempty(h),
   a.h = h;
end;

ann = nneg(a); an = neg(a);
[bnn,bnnlog] = integral(ann,var);
[bn,bnlog] = integral(an,var);
b = tsp(bn+bnn);
if nargout==2,
   blog = pol(zeros(a.s));
end;
if ~isempty(a),
   if nargout<2,
      if any(any(bnlog~=0)) | any(any(bnnlog~=0)),
         warning('Integral is not polynomial, contains logarithm.');
      end;
   else
      if any(any(bnlog~=0)),
         blog = bnlog;
      elseif any(any(bnnlog~=0)),
         blog = bnnlog;
      end;
   end;   
end;

%end .. tsp/integral
