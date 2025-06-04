function [C,D,E,F,degT] = debe(M,arg2);
%DEBE  Deadbeat controllers of discrete-time linear system
%
% The commmand
%   C = DEBE(M)
% computes the polynomial matrix fraction description C of a deadbeat
% controller of the linear discrete-time strictly proper 
% system given by the polynomial matrix fraction M. 
%
% The commmand
%    [NC,DC,E,F,DEGT] = DEBE(M), 
% returns either right
%     C = (NC+E*T)/(DC-F*T) 
% or left
%	  C = (DC-T*F)\(NC+T*E),
% polynomial matrix fraction parametrization of all
% deadbeat controllers depending on the internal
% representation of the matrix fraction M (left/right fraction).
% If the variable of M is 'z' or 'q' then T is a free polynomial
% matrix with entry degrees less than or equal to DEGT. 
% For the variables 'd' and 'z^-1', T is a free polynomial matrix
% such that X{0} is nonsingular.
%
% Note that a negative feedback loop is always considered.
%
% A tolerance TOL can be specified as an additional input argument.
%  
% See also PPLACE, STAB.

%	Author(s): M. Hromcik, M. Sebek 14-10-98
%	Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 1.0 $  $Date: 14-Oct-1998 10:28:34   $
%                         $Date: 04-Jul-2001  J.Jezek   $

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
  [Nc,Dc]=debe(Nm,Dm,opt,tol);
  if opt=='l'
    C=Nc/Dc;
  else
    C=Dc\Nc;
  end;
elseif nargout==5
  [C,D,E,F,degT]=debe(Nm,Dm,opt,tol);
else
  error('Invalid number of outputs. (expecting 1 or 5)');
end;

% this function is just a FRAC wrapper for POL/DEBE