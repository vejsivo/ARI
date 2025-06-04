function [N,D] = testdnd(N,D,hand)
%TESTDND   Test dimensions of polynomials to be
%          a numerator and a denominator of a fraction
%
% The command   [N,D] = TESTDND(N,D,'r')  tests
% dimensions of polynomial matrices N,D, whether they
% be a numerator and a denominator of fraction N*D^-1 .
% Admissible is also D scalar; in this case D is
% multiplied by an eye matrix of coresponding dimension.
%
% The command   [N,D] = TESTDND(N,D,'l')  works
% similarly for fraction  D^-1*N . The default is 'r'.
%
% If not O.K. then an error message is issued.

%     Author: J. Jezek, 07-Aug-2001
%     Copyright(c) 2001 by Polyx, Ltd.

if nargin<3, hand = 'r';
end;

if any(any(isnan(N))) | any(any(isnan(D))) | ...
   any(any(isinf(N))) | any(any(isinf(D))),
      error('Polynomial is not finite.');
end

[Ns1,Ns2] = size(N); [Ds1,Ds2] = size(D);
if Ds1 ~= Ds2,
   error('Denominator matrix is not square.');
end;

if strcmp(hand,'l'),
   if Ds2 ~= Ns1,
      if Ds2 == 1,
         D = diag(repmat(D,Ns1,1));
      else error('Matrices of inconsistent dimensions.');
      end;
   end;
else
   if Ds1 ~= Ns2,
      if Ds1 == 1,
         D = diag(repmat(D,Ns2,1));
      else error('Matrices of inconsistent dimensions.');
      end;
   end;
end;

%end .. private/testdnd
