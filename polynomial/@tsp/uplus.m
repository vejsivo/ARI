function A = uplus(T);
%UPLUS  (+) Unary plus of two-sided polynomial
%              A = +T
%
% See also TSP/UMINUS, TSP/PLUS.

%       Author: J.Jezek 11-8-99
%       Copyright (c) 1999 by Polyx, Ltd.
%
% Effect on other properties:
% A.u: UserData are deleted.
% A.version: set 3.0

A = T;
A.u = [];
A.version = 3.0;

%end .. @tsp/uplus
