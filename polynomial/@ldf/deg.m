function [D,L] = deg(A,varargin)
%DEG     Degree of left-den fraction
%
% Fraction A(v) is treated as a function of complex
% variable 'v' and its Laurent series is taken.
% The degree  DEG(A) shows what is the highest power
% of 'v' in the series. When 'v' is 'z','q','s' or 'p',
% the series is taken in point Inf, when 'z^-1' or 'd',
% in point 0. 
%
% In the command  D = DEG(A,STRING) , the second input
% argument STRING means:
%   D = DEG(A)        default, the same as DEG(A,'mat').
%   D = DEG(A,'mat')  the degree of matrix A.
%   D = DEG(A,'ent')  the matrix of degrees of the A entries.
%   D = DEG(A,'row')  the column vector of row degrees of A.
%   D = DEG(A,'col')  the row vector of column degrees of A.
%
% For fractions in 'z' or 'z^-1', it is possible to specify
% (by an optional input argument VAR) that the degree is to be
% understood by the highest power of 'z' or 'z^-1'.
%
% In the command  [D,L] = DEG(A,STRING) , the second
% output argument L returns the correspoding matrix
% of leading coefficients, the same as first output
% argument in function LCOEF. 
%
% See also POL/DEG, POL/LCOEF, LDF/LCOEF.

%       Author:  J. Jezek, 06-Nov-2002
%       Copyright(c) 2002 by Polyx, Ltd.
%       $ Revision $  $ Date 17-Nov-2002 $

if ~isa(A,'ldf'),
   error('Some argument but not 1st is invalidly ldf.');
end;

li = length(varargin);
if li==0,
   if nargout<2,
      eval('D = deg(mdf(A));','error(peel(lasterr));');
   else
      eval('[D,L] = deg(mdf(A));','error(peel(lasterr));');
   end;
else
   if nargout<2,
      eval('D = deg(mdf(A),varargin{1:li});','error(peel(lasterr));');
   else
      eval('[D,L] = deg(mdf(A),varargin{1:li});','error(peel(lasterr));');
   end;
end;

%end .. @ldf/deg

   
   
   
   
   