function [C,D,E,F,degT] = pplace(M,POLES,arg3);
%PPLACE  Polynomial pole placement
%
% Given the polynomial matrix fraction M and a vector of desired pole
% locations POLES the command
%    C = PPLACE(M,POLES)
% computes the polynomial matrix fraction representation C of a feedback 
% controller for the system M that assigns the closed-loop 
% system poles to the location POLES. The multiplicity of the poles is
% increased if necessary. The resulting system may have real or complex 
% coefficients depending on whether or not the desired poles are self-conjugate.
%
% The commmand
%    [NC,DC,E,F,DEGT] = PPLACE(M,POLES)
% returns either right
%	  C = (NC+E*T)/(DC-F*T)
% or left
%	  C = (DC-T*F)\(NC+T*E)
% polynomial matrix fraction parametrization of all controllers that yield
% the same closed-loop poles depending on the internal representation
% of the matrix fraction M (left/right fraction).
% NC, DC, E and F are polynomial matrices and T is an arbitrary polynomial
% matrix parameter of compatible size and of degree bounded by degT. 
%
% Note that a negative feedback loop is always considered.
%
% NOTE: The user should be aware that for multi-input multi-output (MIMO) 
% systems just assigning pole locations need not be enough. In the MIMO case, 
% the desired behavior typically depends on the closed-loop invariant polynomials 
% rather than on the pole locations only. To assign the desired invariant 
% polynomials, put them into a diagonal matrix R of the same size as D and call 
%    C = pplace(M,R)
% or
%    [NC,DC,E,F,DEGT] = PPLACE(M,R).
%
% A tolerance TOL may be specified as an additional input argument.  

%	Author(s): M. Hromcik, M. Sebek 14-10-98
%	Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 1.0 $  $Date: 14-Oct-1998 10:28:34   $

global PGLOBAL;
eval('PGLOBAL.ZEROING;', 'painit;');
tol = PGLOBAL.ZEROING;

switch nargin
  case {0,1},
    error('Not enough input arguments.');
  
  case 2;
    
  case 3,
    if isnumeric(arg3)
     tol = arg3;
    else 
     error('Invalid 3th argument.');
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
  [Nc,Dc]=pplace(Nm,Dm,POLES,opt,tol);
  if opt=='l'
    C=Nc/Dc;
  else
    C=Dc\Nc;
  end;
elseif nargout==5
  [C,D,E,F,degT]=pplace(Nm,Dm,POLES,opt,tol);
else
  error('Invalid number of outputs. (expecting 1 or 5)');
end;

% this function is just a FRAC wrapper for POL/PPLACE