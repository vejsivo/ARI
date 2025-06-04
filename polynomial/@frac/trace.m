function Rt = trace(R, tol);
%TRACE  Sum of the diagonal elements of fraction
%   
% TRACE(R) is the sum of the diagonal elements of fraction R
% (scalar=den fraction or matrix-den fraction) with zeroing
% activated through the global variable PGLOBAL.ZEROING.
%
% TRACE(R,TOL) works with zeroing specified by the input 
% tolerance TOL. 

%      Author:  J. Jezek  28-Jan-2000
%      Copyright(c) 2000 by Polyx, Ltd.
%      $ Revision $  $ Date 25-Jul-2002 $

global PGLOBAL;

if nargin == 1 | isempty(tol), 
   tol = PGLOBAL.ZEROING;
else
   if ~isa(tol,'double') | length(tol)~=1 | ...
         ~isreal(tol) | tol<0 | tol>1,
      error('Invalid tolerance.');
   end;
end;

Rcl = class(R);
if strcmp(Rcl,'frac'),
   error('Function ''trace'' not defined for variables of class ''frac''.');
end;

Rt = sum(diag(R),[],tol);

%end .. @frac/trace

  