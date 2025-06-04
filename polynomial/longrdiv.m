function [Q,R] = longrdiv(N, D, degree, tol);
%LONGRDIV   Long right divide polynomials
%
% If NUM and DEN are polynomial matrices of consistent dimensions
% in the variables 'z', 'q' or '', then the command
%    [Q,R] = LONGRDIV(NUM,DEN,DEG)
% returns the polynomial matrices 
%    Q(z) = Q1*z + Q2*z^2+ ... + Qn*z^n 
% and 
%	  R(z^-1) = R0 + R1*z^-1 ... RDEG*z^-DEG
% which are the first terms of the Laurent series 
%    L(z) = Qn*z^n +..+ Q1*z + R0 + R1*z^-1 + R2*z^-2 + ...
% of the rational matrix NUM(z)/DEN(z) about the point z = Inf.
%
% If NUM and DEN are polynomial matrices of consistent dimensions
% in the variables 's' or 'p' then the command
%    [Q,R] = LONGLDIV(NUM,DEN,DEG)
% returns the polynomial matrix Q,
%	  Q(s) = Q1*s + Q2*s^2+ ... + Qn*s^n 
% and a three-dimensional array R(:,:,1:DEG+1) such that the first terms 
% of the Laurent series of the rational matrix NUM(s)/DEN(s) about the 
% point s = Inf are
%	  Qn*s^n + .. + Q1*s + R0 + R1*s^-1 + .. + RDEG*s^-DEG,
% with Ri = R(:,:,i+1) for i = 0,1...DEG .
%
% If NUM and DEN are polynomial matrices of consistent dimensions in
% the variables 'd' or 'z^-1' then the command
%    [Q,R] = LONGRDIV(NUM,DEN,DEG)
% returns polynomial matrices Q(z) and R(z^-1) which form the first 
% terms of the Laurent series of NUM(z^-1)/DEN(z^-1) about the point 
% d = z^-1 = 0.
%
% The command
%    Q = LONGRDIV(NUM,DEN,DEG) 
% returns the polynomial part of L(var) including the zero-degree term, 
% that is,
%		R0 + Q1*var + Q2*var^2 + ... + Qn*var^n.
%
% A tolerance TOL may be specified as an optional input argument.
%
% If NUM/DEN is the polynomial matrix fraction description of a 
% discrete-time linear system in the forward shift operator 'z' or 'q'
% or in the backward shift operator 'd' or 'z^-1' then its impulse
% response is G(-N) = QN for -N < 0 and G(N) = RN for N >= 0. 
% If the system is causal then Q is zero polynomial matrix.
%
% For the variables 's' and 'p' the array R contains the first
% Markov parameters of the continuous-time system NUM/DEN. 
% R(:,:,i) is the i-th Markov parameter of the system.
%
% See also: RDIV, LONGLDIV.

%	Author(s): M. Hromcik, M. Sebek 15-9-98
%	Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 8-Oct-1998 10:28:34   $
%       $Revision: 3.0 $  $Date: 16-Feb-2000 15:30:00  J.Jezek  $
%                         $Date: 29-May-2000 14:00:00  J.Jezek  $
%                         $Date: 18-Jul-2000 11:30:00  J.Jezek  $
%                         $Date: 25-Sep-2002 comments  J.Jezek  $
%                         $Date: 28-Feb-2003           J.Jezek  $

global PGLOBAL;
eval('PGLOBAL.ZEROING;', 'painit;');

ni = nargin;
if ni<3,
   degree = [];   
end;
if ni<4, 
   tol = PGLOBAL.ZEROING; 
end;

if nargout == 2,
   eval('[Q,R] = longldiv(N.'', D.'', degree, tol);', ...
      'error(peel(lasterr));');
   Q = Q.';
   if isa(R,'double') & ndims(R)==3,
      R = permute(R,[2,1,3]);
   else
      R = R.';
   end;
else
   eval('Q = longldiv(N.'', D.'', degree, tol);', ...
      'error(peel(lasterr));');
   Q = Q.';
end;  

%end .. longrdiv
