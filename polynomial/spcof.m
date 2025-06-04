function [A, J] = SPCOF(B, varargin)
%SPCOF  Polynomial spectral co-factorization
%
% If B is a continuous-time para-Hermitian polynomial matrix,
% that is, B = B', then
%    [A,J] = SPCOF(B)
% solves the polynomial J-spectral co-factorization problem, i.e.
%     B = AJA',
% where A is Hurwitz and J is a signature matrix.
%
% If B is a discrete-time para-Hermitian polynomial matrix, that is,
%     B = BH' + BH
% and B is positive-definite on the unit circle, then
%     A = SPCOF(B)
% solves the polynomial spectral co-factorization problem, i.e.
%     B = AA',  
% where A is Schur.
%
% For more information see the macro SPF.

%    Authors: D. Henrion, August 10, 1998.
%    Copyright 1998 by Polyx, Ltd.
%    Revised Dec, 2013.

if nargin<1,
   error('Not enough input arguments.');
end;

if nargout < 2,
   eval('A = spf(B.'', varargin);', 'error(peel(lasterr));');
   A = A.';
else
   eval('[A,J] = spf(B.'', varargin);', 'error(peel(lasterr));');
   A = A.';
end;
 
%end .. spcof
