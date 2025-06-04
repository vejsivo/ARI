function B = samp(A,T,tau)
%SAMP    Sample continuous time polynomial
%
% For polynomial in variable 's' or 'p', the command
%   B = SAMP(A,T,TAU)  
% where T,TAU are sampling period (scalar) and
% sampling phase (scalar or vector), T,TAU >= 0,
% in nontrivial cases results in error when any 
% TAU is zero, and returns 0 when all TAUs are
% nonzero.
%
% The arguments T,TAU are optional, the defaults are
%  T = 1, TAU = 0 .
%
% This macro exists only for completeness.
% See also RDF/SAMP, LDF/SAMP, MDF/SAMP, SDF/SAMP.

%         Author: J. Jezek  19-Jun-2000
%         Copyright(c) by Polyx, Ltd.
%         $ Revision $  $ Date 04-Oct-2000 $
%                       $ Date 20-Jan-2001 $

ni = nargin;
if ni<1,
   error('Not enough input arguments.');
end;
if ~isa(A,'pol'),
   error('Some argument but not 1st is invalidly polynomial.');
end;

if ni<=1 | isempty(T),
   T = 1;
else
   if isa(T,'double') & length(T)==1 & isreal(T) & T>=0,
   else
      error('Invalid sampling period.');
   end;
end;

if ni<=2 | isempty(tau),
   tau = 0;
else
   if isa(tau,'double') & ndims(tau)==2 & any(size(tau)==1) & ...
         isreal(tau) & all(tau>=0),
   else
      error('Invalid sampling phase.');
   end;
end;

Av = A.v;
if ~isempty(Av) & ~strcmp(Av,'s') & ~strcmp(Av,'p'),
   error('Invalid variable symbol; must be continuous-time.');
end;

if ~isempty(A) & any(any(A~=0)) & any(tau==0),
   error('In time domain, sampling of Dirac impulse.');
end;

BB = pol(zeros(size(A)));
ltau = length(tau);
if ltau==1,
   B = BB;
else
   B = cell(size(tau));
   for k = 1:ltau,
      B{k} = BB;
   end;
end;

%end .. @pol/samp
