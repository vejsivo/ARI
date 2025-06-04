function C = coef(A,arg2,arg3)
%COEF    Coefficient of constant
%
% For constant matrix A, the command  C = COEF(A,K)
% returns A for K=0, otherwise returns zero matrix of 
% coresponding dimensions.
%
% This macro exists only for completeness.
%
% See also: POL/COEF

%      Author:  J. Jezek, 22-Jul-2002
%      Copyrigth(c) 2002 by Polyx, Ltd.

ni = nargin;
if ni<2,
   error('Not enough input arguments.');
elseif ni==2,
   eval('C = coef(pol(A),arg2);','error(peel(lasterr));');
else
   eval('C = coef(pol(A),arg2,arg3);','error(peel(lasterr));');
end;

%end .. coef


      
