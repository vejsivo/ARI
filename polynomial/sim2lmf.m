function [num,den] = sim2lmf(model, X, U);
%SIM2LMF	Simulink-to-LMF description of a dynamic system.
%
% The command 
%	[NUM,DEN] = SIM2LMF('MODEL')
% finds the linear approximation to the dynamical
% system described by the Simulink scheme named 
% 'MODEL' in the form of right (matrix) fraction.
% Zero initial conditions are considered
% implicitely, but can be added in a vector 
% as an extra input.
%
% See also SIM2RMF, POL/SS, SS2LMF.

%       Author(s): M. Hromcik, M. Sebek 16-2-2000
%       Copyright (c) by PolyX, Ltd.
%       $Revision: 1.0 $ 
%       Modified by J. Jezek, Aug 2001, arg checking

if nargin<1
   error('Not enough input arguments.');
end;
if ~ischar(model)
   error('Argument must contain a string.');
end;

switch nargin,
case 1
   [a,b,c,d] = linmod(model);
case 2
   [a,b,c,d] = linmod(model, X);
case 3
   [a,b,c,d] = linmod(model, X, U);
end;

[num,den] = ss2lmf(a,b,c,d);

%end ..sim2lmf
