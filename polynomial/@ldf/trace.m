function Rt = trace(R,tol);
%TRACE  Sum of the diagonal elements of left-den fraction
%   
% TRACE(R) is the sum of the diagonal elements of R
% with zeroing activated through the global variable
% PGLOBAL.ZEROING.
%
% TRACE(R,TOL) works with zeroing specified by the input 
% tolerance TOL. 

%      Author:  J. Jezek  07-Feb-2000
%      Copyright(c) 2000 by Polyx, Ltd.
%      $ Revision $  $ Date 21-Apr-2000 $
%                    $ Date 28-Feb-2003 $

global PGLOBAL;

if nargin == 1 | isempty(tol), 
   tol = PGLOBAL.ZEROING;
else
   if ~isa(tol,'double'),
      error('Invalid tolerance.');
   end;
end;

eval('R = sdf(R,tol);','error(peel(lasterr));');
Rt = ldf(trace(R,tol));  

%end .. @ldf/trace

  