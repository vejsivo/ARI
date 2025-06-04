function [Q,R] = longldiv(N, D, degree, tol);
%LONGLDIV   Long left divide polynomials
%
% If NUM and DEN are polynomial matrices of consistent dimensions
% in the variables 'z', 'q' or '', then the command
%    [Q,R] = LONGLDIV(NUM,DEN,DEG)
% returns the polynomial matrices 
%	  Q(z) = Q1*z + Q2*z^2+ ... + Qn*z^n 
% and 
%	  R(z^-1) = R0 + R1*z^-1 ... RDEG*z^-DEG
% which are the first terms of the Laurent series 
%	  L(z) = Qn*z^n +..+ Q1*z + R0 + R1*z^-1 + R2*z^-2 + ...
% of the rational matrix DEN(z)\NUM(z) about the point z = Inf.
%
% If NUM and DEN are polynomial matrices of consistent dimensions
% in the variables 's' or 'p' then the command
%    [Q,R] = LONGLDIV(NUM,DEN,DEG)
% returns the polynomial matrix
%	  Q(s) = Q1*s + Q2*s^2+ ... + Qn*s^n 
% and the three-dimensional array R(:,:,1:DEG+1) such that the first 
% terms of the Laurent series of the rational matrix DEN(s)\NUM(s) 
% about the point s = Inf are
%	  Qn*s^n + .. + Q1*s + R0 + R1*s^-1 + .. + RDEG*s^-DEG,
% with  R(:,:,i+1) = Ri  for i = 0,1...DEG .
%
% If NUM and DEN are polynomial matrices of consistent dimensions
% in the variables 'd' or 'z^-1' then the command
%    [Q,R] = LONGLDIV(NUM,DEN,DEG)
% returns the polynomial matrices Q(z) and R(z^-1) or R(d) which
% form the  first terms of the Laurent series of
% DEN(z^-1)\NUM(z^-1) about the point d = z^-1 = 0.
%
% The command
%    Q = LONGLDIV(NUM,DEN,DEG) 
% returns the polynomial part of L(var) including the zero-degree term, 
% that is,
%		R0 + Q1*var + Q2*var^2 + ... + Qn*var^n
%
% A tolerance TOL may be specified as an additional input argument.
%
% If DEN\NUM is the polynomial matrix fraction description of a 
% discrete-time linear system in the forward shift operator 'z' or 'q' 
% or in the backward shift operator 'd' or 'z^-1' then its impulse 
% response is G(-N) = QN for -N < 0 and G(N) = RN for N >= 0. If 
% the system is causal then Q is a zero polynomial matrix.
%
% For the variables 's' and 'p' the 3-d array R contains the first
% Markov parameters of the continuous-time system DEN\NUM. 
% R(:,:,i) is the i-th Markov parameter of the system.
%
% See also: LDIV, LONGRDIV.

%	     Author(s): M. Hromcik, M. Sebek 15-9-98
%	     Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 3.0 $  $Date: 16-Oct-1998 10:28:34   $
%                         $Date: 16-Feb-2000 15:30:00  J.Jezek $ 
%                         $Date: 29-May-2000 11:30:00  J.Jezek $
%                         $Date: 04-Jul-2000 13:00:00  J.Jezek $ 
%                         $Date: 18-Jul-2000 13:15:00  J.Jezek $
%                         $Date: 15-Mar-2002           J.Jezek $
%                         $Date: 25-Sep-2002           J.Jezek $
%                         $Date: 28-Feb-2003           J.Jezek $

global PGLOBAL;
eval('PGLOBAL.ZEROING;', 'painit;');

ni = nargin; no = nargout;
if (no==2 & ni<3) | ni<2,
   error('Not enough input arguments.');
end;
if no==2 & isempty(degree),
   error('Degree argument missing.');
end;
if ni>=3 & ~isempty(degree),
   if ~isa(degree,'double') | length(degree)~=1 | ~isfinite(degree) | ...
         ~isreal(degree) | abs(floor(degree))~=degree,
      error('Invalid 3rd argument; must be nonnegative integer.'); 
   end;
end;
if ni<=3, 
   tol = PGLOBAL.ZEROING;  
else
   if ~isa(tol,'double') | length(tol)~=1 | ~isreal(tol) | ...
         tol<0 | tol>1,
      error('Invalid tolerance.');
   end;
end;

eval('N = pol(N); D = pol(D);', 'error(peel(lasterr));');

Nv = N.v; Dv = D.v;
Ns = N.s; Ds = D.s;

if Ds(1)~=Ds(2), error('Denominator matrix is not square.');
end;
if rank(D, tol) < Ds(1), warning('Denominator matrix is singular.');
end;

%if all(Ds == 1),
%  Ds = [Ns(1) Ns(1)];
%end;  

% Check size:
if Ds(2) ~= Ns(1) & all(Ds - [1 1]),
  error('Matrices of inconsistent dimensions.');
end;

% Variable and sampling period consistency:
[tv,var,N,D] = testvp(N,D);
if tv==0,
   warning('Inconsistent variables.');
elseif tv==2,
   error('Inconsistent variables.');
end;
[th,per,N,D] = testhp(N,D,var);
if ~th,
   warning('Inconsistent sampling periods.');
end;

% d-z conversion if necessary:
isd = strcmp(var,'d') | strcmp(var, 'z^-1');
if isd, [N,D] = reverse(N,D,tol);  % N and D in 'z'
end;

% Make D row reduced if necessary:
if rank(lcoef(D, 'row')) < min(Ds),
  %error('D is not reduced.');
  [D,pom,U] = rowred(D);
  N = U*N;
end;  

% Non-causal part Q - polynomial in 'z':
% eval('[Q, N] = ldiv(N,D, tol);', 'error(lasterr)');
[Q, N] = ldiv(N,D, tol);

Nd = N.d; Dd = D.d;

if no == 2,
	% DEGREE is zero - extract just the quotient:
	if degree == 0,
 	  R = Q{0};
 	  if Q.d >= 1,
 	   Q{0} = 0;
 	  else	% Q.d <= 0,
 	   Q = pol(zeros(Q.s));
 	  end;
 	  return;    
	end;
  
 	% Causal part - long division, polynomial in 'z^-1':
	vr = pol([0 1],1,var);
	shft = degree;	% (degree+D.d-N.d)-1;
   % ldiv(N, D)
	N = shift(N, shft);
   N.v = var; N.h = per;
   [R,pom]=ldiv(N,D,tol);
	%eval('[R, pom] = ldiv(N, D, tol)','R = pol(NaN)*ones(Ns);');%%%%%%
   %Q.v = 'z';
   if ~strcmp(Q.v,'q'), Q.v = 'z';
   end;
	Rc = flipdim(R.c, 3);
	Rc = cat(3, zeros([R.s, shft+1 - size(R.c,3)]), Rc);
	R = pol( Rc(:,:), shft );
   %	R.v = 'z^-1'; R.h = per;
   R.v = 'z^-1';
   if strcmp(var,'s') | strcmp(var,'p'),
      R.h = NaN;
   else
      R.h = per;
   end;
	R = R+Q{0};
  
	if Q.d >= 1,
	  Q{0} = 0;
	else	% Q.d <= 0,
	  Q = pol(zeros(Q.s));
   end;    
  
   if strcmp(var,'d'), R.v = 'd';
   end;
  
	% Continuous case - return the Markov parameters structure:
	if strcmp(var,'s') | strcmp(var,'p'),
	  Q.v = var; 
	  Rs = R.s;
	  R = R.c;
	  
	  %R = permute(R.c, [3,2,1]);
	  %eval('R = R(1:degree,:,:);', ''); 
	  
	  %eval('R = R(:,:,1:degree);', ''); 
     if isempty(R), R = zeros(Rs);
     end;
     degree1 = degree+1;
	  R(:,:,size(R,3)+1:degree1) = zeros([Rs degree1-size(R,3)]);
	end;	
end; % .. if no	

%end .. longldiv
