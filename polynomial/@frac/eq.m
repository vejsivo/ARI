 function c = eq(P,Q,tol);
%EQ (==)   Test if fractions equal
%	    P == Q
%
% The commmand
%    C = (P==Q) 
% performs elementwise comparison between the fractions
% P and Q with a tolerance activated through the global variable
% PGLOBAL.ZEROING. P and Q must have the same dimensions unless
% one is a scalar fraction; in such a case the scalar is compared
% with every entry of the other matrix.
%
% A difference between the variable symbols of P and Q causes a 
% warning message, but does not affect the result.    
% However, if one of the variables is 'z' and the other 'z^-1',
% then the variable names play a role in the comparison;
% no warning is issued in such a case.
% 
% The commmands
%    C = EQ(P,Q) 
% works alike. The commmand
%    C = EQ(P,Q,TOL) 
% works with tolerance specified by the input tolerance TOL.
%  
% See also: FRAC/NE

%      Author:  J. Jezek, 19-Sep-2001
%      Copyright(c) 2001 by Polyx, Ltd.
%      $ Revision $   $ Date 29-Jan-2002 $

global PGLOBAL;

ni = nargin;
if ni<=1,
   error('Not enough input arguments.');
elseif ni==2 | isempty(tol),
   tol = PGLOBAL.ZEROING;
else
   if ~isa(tol,'double'),
      error('Invalid tolerance.');
   end;
end;

eval('P = mdf(P);', 'error(''Invalid 1st argument.'');');
eval('Q = mdf(Q);', 'error(''Invalid 2nd argument.'');');

[tv,Rv,P,Q] = testvf(P,Q);
if ~tv, warning('Inconsistent variables.');
end;
[th,Rh,P,Q] = testhf(P,Q,Rv);
if ~th, warning('Inconsistent sampling periods.');
end;

eval('c = eq(times(P.num,Q.den,tol),times(P.den,Q.num,tol),tol);',...
   'error(peel(lasterr));');

%end .. @frac/eq
