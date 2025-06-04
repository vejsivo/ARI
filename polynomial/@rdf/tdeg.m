function [D,L] = tdeg(A,varargin)
%TDEG    Trailing degree of right-den fraction
%
% Fraction A(v) is treated as a function of complex
% variable 'v' and its Laurent series is taken.
% The trailing degree TDEG(A) shows what is the
% lowest power of 'v' in the series. When 'v' is 'z','q',
% 's' or 'p', the series is taken in point Inf,
% when 'z^-1' or 'd', in point 0.
%
% In the command  D = TDEG(A,STRING) , the second input
% argument STRING means:
%   D = TDEG(A)        default, the same as TDEG(A,'mat').
%   D = TDEG(A,'mat')  the trailing degree of matrix A.
%   D = TDEG(A,'ent')  the matrix of trailing degrees of the A entries.
%   D = TDEG(A,'row')  the column vector of row trailing degrees of A.
%   D = TDEG(A,'col')  the row vector of column trailing degrees of A.
%
% For fractions in 'z' or 'z^-1', it is possible to specify
% (by an optional input argument VAR) that the triling degree is to be
% understood by the highest power of 'z' or 'z^-1'.
%
% In the command  [D,L] = TDEG(A,STRING) , the second
% output argument L returns the correspoding matrix
% of trailing coefficients, the same as first output
% argument in function TCOEF.
%
% See also POL/TDEG, POL/TCOEF, RDF/TCOEF.

%       Author:  J. Jezek, 07-Nov-2002
%       Copyright(c) 2002 by Polyx, Ltd.
%       $ Revision $  $ Date 17-Nov-2002 $

if ~isa(A,'rdf'),
   error('Some argument but not 1st is invalidly rdf.');
end;

li = length(varargin);
if li==0,
   if nargout<2,
      eval('D = tdeg(mdf(A));','error(peel(lasterr));');
   else
      eval('[D,L] = tdeg(mdf(A));','error(peel(lasterr));');
   end;
else
   if nargout<2,
      eval('D = tdeg(mdf(A),varargin{1:li});','error(peel(lasterr));');
   else
      eval('[D,L] = tdeg(mdf(A),varargin{1:li});','error(peel(lasterr));');
   end;
end;

%end .. @rdf/tdeg

