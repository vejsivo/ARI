function c = ne(P,Q,tol);
%NE (~=)   Test if fractions unequal
%	    P ~= Q
%
% The commmand
%    C = (P~=Q) 
% performs elementwise comparison between the fractions P and Q,
% using a tolerance activated through PGLOBAL.ZEROING.
% P and Q must have the same dimensions unless one is a scalar;
% in such a case the scalar is compared with every element
% of the other matrix.
%
% A difference between the variable symbols of P and Q causes
% a warning message, but does not affect the result.    
% However, if one of the variables is 'z' and the other 'z^-1',
% then the variable names play a role in the comparison;
% no warning is issued in such a case.
% 
% The commmands
%    C = NE(P,Q) 
% works alike. The commmand
%    C = NE(P,Q,TOL) 
% works with tolerance specified by the input tolerance TOL.
%  
% See also: FRAC/EQ.

%       Author:  J. Jezek  26-Jan-2000
%       Copyright(c) 2000 by Polyx, Ltd.
%       Modified  19-Sep-2001

global PGLOBAL;

ni = nargin;
if ni<=1,
   error('Not enough input arguments.');
elseif ni==2 | isempty(tol),
   tol = PGLOBAL.ZEROING;
elseif ~isa(tol,'double'),
   error('Invalid tolerance.');
end;

eval('c = ~eq(P,Q,tol);','error(peel(lasterr));');

%end .. @frac/ne
