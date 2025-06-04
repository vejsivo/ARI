function [L,D] = tcoef(A,varargin)
%TCOEF   Trailing coefficient of right-den fraction
%
% Fraction A(v) is treated as a function of complex
% variable 'v' and its Laurent series is taken.
% The trailing coefficient TCOEF(A) is the coefficient
% with the lowest power of 'v' in the series.
% When 'v' is 'z','q','s' or 'p', the series is taken
% in point Inf, when 'z^-1' or 'd', in point 0.
%
% In the command  D = TCOEF(A,STRING) , the second input
% argument STRING means:
%   D = TCOEF(A)        default, the same as TCOEF(A,'mat').
%   D = TCOEF(A,'mat')  the trailing coefficient matrix of A.
%   D = TCOEF(A,'ent')  the matrix of scalar trailing coefficients
%                         of A entries.
%   D = TCOEF(A,'row')  the row trailing coefficient matrix of A.
%   D = TCOEF(A,'col')  the column trailing coefficient matrix of A.
%
% For fractions in 'z' or 'z^-1', it is possible to specify
% (by an optional input argument VAR) that the trailing coefficient
% is to be understood by the lowest power of 'z' or 'z^-1'.
%
% In the command  [L,D] - TCOEF(A,STRING) , the second
% output argument D returns the corresponding matrix or vector of
% trailing degrees, the same as the first output argument in
% function TDEG. In the case if infinite series, the
% corresponding degree is -Inf.
%
% See also POL/TCOEF, POL/TDEG, RDF/TDEG.

%       Author:  J. Jezek, 07-Nov-2002
%       Copyright(c) 2002 by Polyx, Ltd.
%       $ Revision $  $ Date 17-Nov-2002 $

if ~isa(A,'rdf'),
   error('Some argument but not 1st is invalidly rdf.');
end;

li = length(varargin);
if li==0,
   if nargout<2,
      eval('L = tcoef(mdf(A));','error(peel(lasterr));');
   else
      eval('[L,D] = tcoef(mdf(A));','error(peel(lasterr));');
   end;
else
   if nargout<2,
      eval('L = tcoef(mdf(A),varargin{1:li});','error(peel(lasterr));');
   else
      eval('[L,D] = tcoef(mdf(A),varargin{1:li});','error(peel(lasterr));');
   end;
end;

%end .. @rdf/tcoef

   