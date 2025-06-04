function [Q,R] = imag(P,tol)
%IMAG    Imaginary part of fraction
%
% For fraction P (scalar-den fraction or matrix-den fraction), the
% command  Q = IMAG(P)  returns fraction Q, the imaginary part of P.
% As the computation is not trivial (it includes addition and
% multiplication of polynomials), the zeroing tolerance can be
% specified by an optional input argument TOL. In the optional
% output argument R, the real part can be also returned.
%
% See also RDF/REAL, LDF/REAL, FRAC/REAL, RDF/ISREAL, LDF/ISREAL,
%          MDF/ISREAL, SDF/ISREAL.

%        Author:  J. Jezek  26-Jan-2000
%        Copyright(c) 2000 by Polyx, Ltd.
%        Revision  $ Date: 21-Sep-2001 $
%                  $ Date: 25-Jul-2002 $
%                  $ Date: 30-Sep-2002 $
%                  $ Date: 14-Oct-2002 $

global PGLOBAL;

if nargin==1,
   tol = PGLOBAL.ZEROING;
else
   if ~isa(tol,'double') | (~isempty(tol) & ...
         (length(tol)~=1 | ~isreal(tol) | tol<0 | tol>1)),
      error('Invalid tolerance.');
   end;
end;

Pcl = class(P);
if strcmp(Pcl,'frac'),
   error('Function ''imag'' not defined for variables of class ''frac''.');
end;

if isreal(P.num) & isreal(P.den),
   Q = P.*0;
else
   A = real(P.num); B = imag(P.num);
   C = real(P.den); D = imag(P.den);
   den = plus(times(C,C,tol),times(D,D,tol),tol);
   Qn = minus(times(B,C,tol),times(A,D,tol),tol);
   Q = P; Q.num = Qn; Q.den = den;
   if strcmp(P.p,'prop'), Qp = 'prop'; Qtp = P.tp;
   else Qp = 'prop?'; Qtp = [];
   end;
   Qc = 'cop?'; Qtc = []; Qr = 'red?';
   props(Q,Qp,Qtp,Qc,Qtc,Qr);
end;

if strcmp(PGLOBAL.COPRIME,'cop'),
   Q = coprime(Q,tol);
end;
if strcmp(PGLOBAL.REDUCE,'red'),
   Q = reduce(Q,tol);
else
   Q = smreduce(Q);
end;
if strcmp(PGLOBAL.DEFRACT,'defr'),
   Q = defract(Q);
end;

if nargout==2,
   if isreal(P.num) & isreal(P.den),
      R = P;
   else
      Rn = plus(times(A,C,tol),times(B,D,tol),tol);
      R = P; R.num = Rn; R.den = den;
      if strcmp(P.p,'prop'), Rp = 'prop'; Rtp = P.tp;
      else Rp = 'prop?'; Rtp = [];
      end;
      Rc = 'cop?'; Rtc = []; Rr = 'red?';
      props(R,Rp,Rtp,Rc,Rtc,Rr);
   end;
   
   if strcmp(PGLOBAL.COPRIME,'cop'),
      R = coprime(R,tol);
   end;
   if strcmp(PGLOBAL.REDUCE,'red'),
      R = reduce(R,tol);
   else
      R = smreduce(R);
   end;
   if strcmp(PGLOBAL.DEFRACT,'defr'),
      R = defract(R);
   end;
end;

%end .. @frac/imag
