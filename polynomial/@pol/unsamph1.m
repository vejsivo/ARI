function  F = unsamph1(G,T,varargin)
%UNSAMPH1     Un- sample and first order hold discrete
%                           time polynomial
%
% Polynomial G is supposed to be a result of sampling
% and first order holding of continuous time polynomial
% F with period T and phase TAU. The command
%     F = UNSAMPH1(G,T,TAU)
% restores the original polynomial F.
%
% For polynomial G in variable z^-1 or d, its degree 
% must not be greater than 1. In variable z or q,
% only constant polynomial is admissible. The degree
% of resulting polynomial F is not greater than 1.
% This holds for all T and TAU. All other cases result
% in error. Argument T is compulsory, TAU is optional,
% default 0.
%
% This macro exists only for completeness.
% See also FRAC/UNSAMPH1.

if nargin<2,
   error('Not enough input arguments.');
end;
if ~isa(G,'pol'),
   error('Invalid 1st argument.');
end;

lv = length(varargin);
eval('F = unsamph1(sdf(G),T,varargin{1:lv});', ...
   'error(peel(lasterr));');
F = defract(F);

%end .. @pol/unsamph1
