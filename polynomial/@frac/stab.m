function [C,D,E,F] = stab(M,arg2);
%STAB  Stabilizing controllers for a linear system
%
% The commmand 
%    C = STAB(M) 
% computes the polynomial matrix fraction description C
% of a stabilizing controller for the linear system given by
% a coprime polynomial matrix fraction M. 
%
% The commmand
%    [NC,DC,E,F] = STAB(M)
% returns either right
%     C = Y/X = (NC*P+E*T)/(DC*P-F*T)
% or left
%	  C = X\Y = (P*DC-T*F)\(P*NC+T*E)
% polynomial matrix fraction parametrization of all
% stabilizing controllers depending on the internal
% representation of the matrix fraction M (left/right fraction).
% NC, DC, E and F are polynomial matrices, T is an arbitrary
% polynomial matrix parameter and P is a stable polynomial 
% matrix such that Y/X (or X\Y, respectively) is proper.
%
% Note that a negative feedback loop is always considered.
%
% A tolerance TOL may be specified as an additional input argument.
%  
% See also PPLACE, DEBE.

%	Author(s): M. Hromcik, M. Sebek 14-10-98
%	Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 1.0 $  $Date: 14-Oct-1998 10:28:34   $

global PGLOBAL;
eval('PGLOBAL.ZEROING;', 'painit;');
tol = PGLOBAL.ZEROING;

switch nargin
   
  case 0,
    error('Not enough input arguments.');
    
  case 1,
    
  case 2,
    if isnumeric(arg2) 
     tol = arg2;
    else 
     error('Invalid 2nd argument.');
  end;
  
  otherwise error('Too many input arguments.');
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
  [Nc,Dc]=stab(Nm,Dm,opt,tol);
  if opt=='l'
    C=Nc/Dc;
  else
    C=Dc\Nc;
  end;
elseif nargout==4
  [C,D,E,F]=stab(Nm,Dm,opt,tol);
else
  error('Invalid number of outputs. (expecting 1 or 4)');
end;

% this function is just a FRAC wrapper for POL/STAB