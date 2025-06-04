function Pshift = shift(P,n,var);
%SHIFT  Shift a constant
%
% The commmand
%    PS = SHIFT(P,N) 
% with constant (i.e. standard Matlab) marix P and
% with integer N, computes  PS(VAR) = P * VAR^N,
% where VAR is either the standard polynomial variable
% or it is given by an optional input argument VAR.
%
% This macro exists only for completeness,
% see also POL/SHIFT.

%	Author(s): M. Hromcik, M. Sebek 20-5-98
%	Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 10-Jun-1998 10:28:34   $
%                         $Date: 13-Aug-2001  J.Jezek   $
global PGLOBAL;
eval('PGLOBAL.ZEROING;', 'painit;'); 

na = nargin;
if na<2,
   error('Not enough input arguments.');
end;
if ~isa(n,'double'),
   error('Invalid shift; must be integer.');
end;

if na==2,
   eval('Pshift = shift(pol(P),n);', 'error(peel(lasterr));');
else
   if ~(isa(var,'char') | isa(var,'pol')),
      error('Invalid variable symbol.');
   end;
   eval('Pshift = shift(pol(P),n,var);', 'error(peel(lasterr));');
end;

%end .. shift
