function C = plqg(M, Q1,R1, Q2,R2, arg6);
%PLQG  Linear-quadratic-Gaussian (LQG) controller design using polynomial methods
%	
% The commmands
%    C = PLQG(M,Q1,R1,Q2,R2)
% compute the LQC-optimal dynamic output feedback controller C 
% for the plant with proper transfer matrix M (see Fig. 1). 
% Q1 and R1 are the covariances of the input and output white noises. 
% The controller minimizes the steady-state value of 
%     E{y'*Q2*y + u'*R2* u}
%
%  	Fig. 1: The closed loop system.
%
%	       Input noise         Output noise
%	     (covariance Q1)      (covariance R1)
%	            |                    |
%	        u   |     +-------+      |   y 
% 	      +---->o-->--|   M   |--->--o-----+
%	      |           +-------+            |
% 	      |   |<------- PLANT ------->|    |
%	      |                                |
% 	      |         +------------+         |
%	      +----<----|    - C     |----<----+
% 	                +------------+   
%	             |<-- CONTROLLER -->|
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.
%
% NOTE: If the plant is continuous time and is not strictly proper
% then the internal properness of the closed loop system is not guaranteed.
%
% See also DEBE, STAB, PPLACE.

%       Author(s): M. Hromcik, M. Sebek 11-11-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 1.0 $  $Date: 11-Nov-1998 10:28:34   $
%                         $Date: 25-Aug-2001  J.Jezek  arg check, sample per $

global PGLOBAL;
eval('PGLOBAL.ZEROING;', 'painit;');
tol = PGLOBAL.ZEROING;

switch nargin
  case {0,1,2,3,4},
    error('Not enough input arguments.');
  
  case 5;
    
  case 6,
    if isnumeric(arg6)
     tol = arg3;
    else 
     error('Invalid 6th argument.');
    end;
    
  otherwise
    error('Too many input arguments.');
end;	% switch  

if length(tol)~=1 | ~isreal(tol) | tol<0 | tol>1,
   error('Invalid tolerance.');
end;

opt = '-';

if isa(M,'ldf')
  opt='l';
elseif isa(M,'rdf')
  opt='r';
else
  error('Invalid 1st argument. (expecting LDF or RDF)');
end;

Nm=M.num;
Dm=M.den;

if nargout<=1
  [Nc,Dc]=plqg(Nm,Dm,Q1,R1,Q2,R2,opt,tol);
  if opt=='l'
    C=Nc/Dc;
  else
    C=Dc\Nc;
  end;
else
  error('Too many outputs.');
end;

% this function is just a FRAC wrapper for POL/PLQG