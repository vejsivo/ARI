function  F = unsamph1(G,T,varargin)
%UNSAMPH1     Un- sample and first order hold discrete
%                            time constant
%
% Constant G (standard MATLAB matrix) is supposed to be
% a result of sampling and first order holding of
% continuous time polynomial F with period T and
% phase TAU. The command   F = UNSAMPH1(G,T,TAU)
% restores the original polynomial F.
%
% The degree of F is not greater than 1, for all T
% and TAU. Argument T is compulsory, TAU is optional,
% default 0.
%
% This macro exists only for completeness.
% See also  FRAC/UNSAMPH1.

%       Author: J.Jezek  01-Dec-2000
%       Copyright(c) 2000 by Polyx, Ltd.

if nargin<2,
   error('Not enough input arguments.');
end;
if ~isa(G,'double'),
   error('Invalid first argument.');
end;

lv = length(varargin);
eval('F = unsamph1(sdf(G),T,varargin{1:lv});', ...
   'error(peel(lasterr));');
F = declass(F);

%end .. unsamph1
