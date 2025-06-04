function [Q,C] = declass(R)
%DECLASS     Declass fraction (convert to possible lower class)
%
% The command  Q = DECLASS(R)  declasses fraction R, i.e. converts
% it to two-sided polynomial (only in 'z' or 'z^-1'), polynomial
% or double (standard Matlab matrix) if possible.
%
% The command  [Q,C] = DECLASS(R)
% returns also the resulting class in C.

% See also FRAC/TSP, FRAC/POL, FRAC/DOUBLE.

%      Author:  J. Jezek  28-Jan-2000
%      Copyright(c) 2000 by Polyx, Ltd.
%      $ Revision $  $ Date 25-Jul-2002 $
%                    $ Date 14-Oct-2002 $
%                    $ Date 07-Nov-2002  bug $
%                    $ Date 12-Nov-2002 $

cl = class(R);
if strcmp(cl,'frac'),
   error('Function ''declass'' not defined for variables of class ''frac''.');
end;

%end .. @frac/declass

