function [L,D] = lcoef(A,varargin)
%LCOEF   Leading coefficient of right-den fration
%
% Fraction A(v) is treated as a function of complex
% variable 'v' and its Laurent series is taken.
% The leading coefficient LCOEF(A) is the coefficient
% with the highest power of 'v' in the series.
% When 'v' is 'z','q','s' or 'p', the series is taken
% in point Inf, when 'z^-1' or 'd', in point 0.
%
% In the command  L = LCOEF(A,STRING) , the second
% input argument STRING means:
%   L = LCOEF(A)        default, the same as LCOEF(A,'mat').
%   L = LCOEF(A,'mat')  the leading coefficient matrix of A.
%   L = LCOEF(A,'ent')  the matrix of scalar leading coefficients
%                         of A entries.
%   L = LCOEF(A,'row')  the row leading coefficient matrix of A.
%   L = LCOEF(A,'col')  the column leading coefficient matrix of A.
%
% For fractions in 'z' or 'z^-1', it is possible to specify
% (by an optional input argument VAR) that the leading coefficient
% is to be understood by the highest power of 'z' or 'z^-1'.
%
% In the command  [L,D] - LCOEF(A,STRING) , the second
% output argument D returns the corresponding matrix or vector of
% degrees, the same as the first output argument in function DEG.
%
% See also POL/LCOEF, POL/DEG, RDF/DEG.

%       Author:  J. Jezek, 07-Nov-2002
%       Copyright(c) 2002 by Polyx, Ltd.
%       $ Revision $  $ Date 17-Nov-2002 $

if ~isa(A,'rdf'),
   error('Some argument but not 1st is invalidly rdf.');
end;

li = length(varargin);
if li==0,
   if nargout<2,
      eval('L = lcoef(mdf(A));','error(peel(lasterr));');
   else
      eval('[L,D] = lcoef(mdf(A));','error(peel(lasterr));');
   end;
else
   if nargout<2,
      eval('L = lcoef(mdf(A),varargin{1:li});','error(peel(lasterr));');
   else
      eval('[L,D] = lcoef(mdf(A),varargin{1:li});','error(peel(lasterr));');
   end;
end;

%end .. @rdf/lcoef
