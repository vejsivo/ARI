function [G,U] = gld(varargin)
%GLD Greatest left divisor
%
% The command
%    G = GLD(N1, N2, .., Nk) 
% computes a greatest left divisor G of several polynomial matrices 
% N1, N2, .., Nk (k > 0) that all have the same numbers of rows.
%
% The columns of G form a polynomial basis for the module spanned by
% the columns of N0 = [N1 N2 .. Nk]. Note that this basis is not
% necessarily of minimum degree. The number of rows in G is equal to
% the rank of N0. If N0 has full row rank then G is square.
%
% The command
%    GLD(N1, .., 'gau') 
% computes the divisor through a modified version of Gaussian elimination. 
% This method is preferable esthetically and generally results in a divisor 
% of low degree. However, it may not be numerically reliable. This is 
% the default method.
%
% The function call
%    GLD(N1, .., 'syl') 
% computes the divisor through stable reductions of Sylvester matrices. 
% This method is preferable numerically but may result in a divisor of
% high degree.
%
% Matrices Mi such that Ni = G*Mi may be recovered with Mi = AXB(G,Ni).
%
% The command
%    [G,U] = GLD(N1, ..) 
% additionally returns a unimodular matrix U such that N0*U = [ G 0 .. 0 ].  
% For two input matrices N1 and N2 U may be split into two pairs of right 
% coprime polynomial matrices (P,Q) and (R,S) such that
%
%    U = [ P R |      N1*P + N2*Q = G
%        | Q S ]      N1*R + N2*S = 0.
%
% A tolerance TOL may be specified as an additional input argument.
%
% See also: GRD.

%     Author: D. Henrion, September 14, 1998.
%     Updated to 3.0 by D. Henrion, August 30, 2000.
%     Copyright 1998-2000 by Polyx, Ltd.

% The function is dual to GRD.

% Transpose input arguments.

if nargin==0,
   error('Not enough input arguments.');
end;
argin = {};
for i = 1:nargin,
 var = varargin{i};
 if isa(var, 'char'), argin{i} = var;
 else argin{i} = var.'; end;
end;

% Call GRD.

if nargout < 2,
 eval('G = grd(argin).'';', 'error(peel(lasterr));');
else
 eval('[G,U] = grd(argin);', 'error(peel(lasterr));');
 U = U.'; G = G.';
end;

%end .. gld
