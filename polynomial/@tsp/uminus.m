function A = uminus(T);
%UMINUS  Unary minus of two-sided polynomial
%              A = -T
%
% See also TSP/UPLUS, TSP/MINUS.

%       Author: J.Jezek, 11-8-99
%       Copyright (c) 1999 by Polyx, Ltd.
%       $ Revision $  $ Date 06-Feb-2002 $

% Effect on other properties:
% A.version: set 3.0

A = T; A.p = -T.p;
A.u = []; A.version = 3.0;

%end .. @tsp/uminus
