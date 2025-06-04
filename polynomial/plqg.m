function [Nc,Dc] = plqg(B,A, Q1,R1, Q2,R2, arg7, arg8);
%PLQG  Linear-quadratic-Gaussian (LQG) controller design using polynomial methods
%	
% The commmands
%    [NC,DC] = PLQG(N,D,Q1,R1,Q2,R2)
%    [NC,DC] = PLQG(N,D,Q1,R1,Q2,R2,'l')
% compute the LQC-optimal dynamic output feedback controller NC/DC 
% for the plant with proper transfer matrix D\N (see Fig. 1). 
% Q1 and R1 are the covariances of the input and output white noises. 
% The controller minimizes the steady-state value of 
%     E{y'*Q2*y + u'*R2* u}
%
%  	Fig. 1: The closed loop system.
%
%	       Input noise         Output noise
%	     (covariance Q1)      (covariance R1)
%	            |                    |
%	        u   |     +-------+      |   y 
% 	      +---->o-->--|  D\N  |--->--o-----+
%	      |           +-------+            |
% 	      |   |<------- PLANT ------->|    |
%	      |                                |
% 	      |         +------------+         |
%	      +----<----|   -NC/DC   |----<----+
% 	                +------------+   
%	             |<-- CONTROLLER -->|
%
% The commmand
%    [NC,DC] = PLQG(N,D, Q1,R1, Q2,R2, 'r') 
% returns such a compensator DC\NC that is optimal for the system 
% given in the form of right polynomial matrix fraction N/D. 
%
% A tolerance TOL may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.
%
% NOTE: If the plant is continuous time and is not strictly proper
% then the internal properness of the closed loop system is not guaranteed.
%
% See also DEBE, STAB, PPLACE.

%       Author(s): M. Hromcik, M. Sebek 11-11-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 1.0 $  $Date: 11-Nov-1998 10:28:34   $
%                         $Date: 25-Aug-2001  J.Jezek  arg check, sample per $

global PGLOBAL;
eval('PGLOBAL.ZEROING;', 'painit;');

tol = PGLOBAL.ZEROING;
opt = 'l';

if nargin<6, error('Not enough input arguments.');
else
  switch nargin
    case 7,
      if isnumeric(arg7), 
         tol = arg7;
      elseif ischar(arg7), 
         opt = arg7;
      else
         error('Invalid 7th argument.');
      end;
    case 8,
      if isnumeric(arg7), 
         tol = arg7;
      elseif ischar(arg7),
         opt = arg7;
      else
         error('Invalid 7th argument.');
      end;
      if isnumeric(arg8), 
         tol = arg8;
      elseif ischar(arg8), 
         opt = arg8;
      else
         error('Invalid 8th argument.');
      end;
  end;            
end;

if ~isempty(tol)
   if length(tol)~=1 | ~isreal(tol) | tol<0 | tol>1,
      error('Invalid tolerance.');
   end;
end;

eval('A = pol(A); B = pol(B);', ...
   'error(peel(lasterr));');

verbose = strcmp(PGLOBAL.VERBOSE,'yes');

% Variables and sampling periods:
[tv,var,B,A] = testvp(B,A);
if tv==2,
   [th,h,B,A] = testhp(B,A,var);
   if ~th,
      warning('Inconsistent sampling periods.');
   end;
   if strcmp(A.v,'z^-1'),
      dg = deg(A);
      B = shift(B,dg); A = rev(A,dg);
      A.v = 'z';
   else
      dg = deg(B);
      B = rev(B,dg); A = shift(A,dg);
      B.v = 'z';
   end;
   var = 'z';
elseif ~tv,
   warning('Inconsistent variables.');
end;

% z^-1 - to - z if needed:
revrsd = 0;
if strcmp(var, 'z^-1') | strcmp(var, 'd'),
   [B,A] = reverse(B,A,tol,opt);
   B.v = 'z'; A.v = 'z';
   revrsd = 1;
elseif strcmp(var,'q'),
   B.v = 'z'; A.v = 'z';
end;

[th,h,B,A] = testhp(B,A,var);
if ~th,
   warning('Inconsistent sampling periods.');
end;

% Sizes:
As = A.s; Bs = B.s;
if As(1)~=As(2), error('Denominator matrix is not square.'); end;

switch opt,

  case 'l',
     if As(1) ~= Bs(1),
        if As(1) == 1,
           A = diag(repmat(A,1,Bs(1)));
           As(1) = Bs(1); As(2) = Bs(1);
        else
           error('Matrices of inconsistent dimensions.');
        end;
    end;  

    isscalar = all(Bs==1);
    A1 = A; B1 = B;

    if isscalar, 
      B2 = B1; A2 = A1;
    else, 
      if rank(lcoef(A1,'row')),
        [A1,aux,U] = rowred(A1);
        B1 = U*B1;
      end;  
      [B2,A2] = lmf2rmf(B1,A1,tol); 
    end;
    
    if ~isproper(B1,A1,'l'),
      error('The plant is not proper.');
    end;  
    
  case 'r',
     if As(2) ~= Bs(2),
        if As(2) == 1,
           A = diag(repmat(A,1,Bs(2)));
           As(1) = Bs(2); As(2) = Bs(2);
        else
           error('Matrices of inconsistent dimensions.');
        end;
    end;   

    isscalar = all(Bs==1);
    A2 = A; B2 = B;

    if isscalar, 
      B1 = B; A1 = A;
    else, 
      if rank(lcoef(A2,'col')),
        [A2,aux,U] = colred(A2);
        B2 = B2*U;
      end;  
      [B1,A1] = rmf2lmf(B2,A2,tol); 
    end;

    if ~isproper(B2,A2,'r'),
      error('The plant is not proper.');
    end;  
    
  otherwise error('Invalid command option.');
  
end;    

% Q and R symmetric, non-negative definite; sizes:
if size(Q1,1) ~= size(Q1,2) | ~isequal(Q1,Q1'),
  error('The covariance matrix of output noise must be square and symmetric.'); 
end;
if size(R1,1) ~= size(R1,2) | ~isequal(R1,R1'),
  error('The covariance matrix of input noise must be square and symmetric.'); 
end;
if size(R2,1) ~= size(R2,2) | ~isequal(R2,R2'),
  error('The weighting matrix R2 must be square and symmetric.'); 
end;
if size(Q2,1) ~= size(Q2,2) | ~isequal(Q2,Q2'),
  error('The weighting matrix Q2 must be square and symmetric.'); 
end;

%if any(real(eig(R1))<=0), error('R1 must be positive definite.'); end;
%if any(real(eig(R2))<=0), error('R2 must be positive definite.'); end;

% Sizes of Q's ,R's:
if size(Q1,1)~=Bs(2), error('Inconsistent dimensions of Q1 and N.'); end;
if size(Q2,1)~=Bs(1), error('Inconsistent dimensions of Q2 and N.'); end;
if size(R1,1)~=Bs(1), error('Inconsistent dimensions of R1 and N.'); end;
if size(R2,1)~=Bs(2), error('Inconsistent dimensions of R2 and N.'); end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For discrete case foot the degrees:
if any(strcmp(var,{'z','q','d','z^-1'})),

  % Spectral factors C1',C2' [cf. 1, p. 306] in z^-1, then C1,C2 in z: 

  %%%%%%%%%%%%%%%%%%%%%%
  sft1 = A1.d - B1.d;
  fac1 = B1*Q1*shift(B1',sft1) + A1*R1*A1'; fac1.v = 'z^-1';
  if verbose, disp('The 1st spectral factorization running ...'); end;
  C1 = spcof(fac1, tol); 	% = C1' in z^-1
  degrowC1 = deg(C1,'row');
  C1 = (C1.')';	C1.v = 'z';	
  C1d = C1.d;
  for i=1:length(degrowC1),
    C1(i,:) = shift(C1(i,:), -C1d+degrowC1(i));
  end;  
  degrowC1 = deg(C1,'row');
  if verbose, disp('The 1st spectral factorization finished.'); end;
  
  sft2 = A2.d - B2.d;
  fac2 = shift(B2',sft2)*Q2*B2 + A2'*R2*A2; fac2.v = 'z^-1';
  if verbose, disp('The 2nd spectral factorization running ...'); end;
  C2 = spf(fac2, tol);
  degcolC2 = deg(C2,'col');
  C2 = (C2.')';	C2.v = 'z';
  C2d = C2.d;
  for i=1:length(degcolC2),
    C2(:,i) = shift(C2(:,i), -C2d+degcolC2(i));
  end;  
  degcolC2 = deg(C2,'col');
  if verbose, disp('The 2nd spectral factorization finished.'); end;
  %%%%%%%%%%%%%%%%%%%%%%%

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % C11 = spcof(B1*Q1*B1','red');
 % C12 = spcof(A1*R1*A1','red');
 % C1 = [C11 C12];
 % 
 % C21 = spf(B2'*Q2*B2);
 % C22 = spf(A2'*R2*A2);
 % C2 = [C21;C22];
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Foot the degrees:
  degrowA1 = deg(A1,'row');
  if any(degrowA1 - degrowC1),
     C1 = diag( pol([0 1],1,'z').^((degrowA1-degrowC1)) ) * C1;
     C1.h = h;
  end  
  
  degcolA2 = deg(A2,'col');
  if any(degcolA2 - degcolC2),
     C2 = C2*diag( pol([0 1],1,'z').^((degcolA2-degcolC2)) );
     C2.h = h;
  end  

else,

  % Spectral factors C1, C2:
  fac1 = [B1 A1]*[Q1 zeros([size(Q1,1) size(R1,2)]);zeros([size(R1,1) size(Q1,2)]) R1]*[B1 A1]';
  if verbose, disp('The 1st spectral factorization running ...'); end;
  C1 = spcof(fac1, tol);
  if verbose, disp('The 1st spectral factorization finished.'); end;

  fac2 = [B2; A2]'*[Q2 zeros([size(Q2,1) size(R2,2)]);zeros([size(R2,1) size(Q2,2)]) R2]*[B2; A2];
  if verbose, disp('The 2nd spectral factorization running ...'); end;
  C2 = spf(fac2, tol);
  if verbose, disp('The 2nd spectral factorization finished.'); end;

  if ~isproper(B2,A2,'r', 'strictly'),
      warning('The plant is not strictly proper. The closed loop system may not be internally proper.');
  end;  

end;
                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isscalar,
  [X2,Y2] = axbyc(A1,B1,C1,tol);
  T2 = ldiv(C2*Y2, A2, tol);  
  Dc = X2*C2 + B2*T2;
  Nc = Y2*C2 - A2*T2;

else,
  switch opt,
    case 'l',
      [X2,Y2] = axbyc(A1,B1,C1,tol);
      [U,V] = rmf2lmf(C2,A2,tol);
      T2 = ldiv(U*Y2, V, tol);
      [T2, C2] = lmf2rmf(T2, C2, tol);
      Dc = X2*C2 + B2*T2;
      Nc = Y2*C2 - A2*T2; 
     
    case 'r',
      [X1,Y1] = xaybc(A2,B2,C2,tol);
      [U,V] = lmf2rmf(C1,A1,tol);
      T1 = rdiv(Y1*U, V, tol);
      [T1, C1] = rmf2lmf(T1, C1, tol);
      Dc = C1*X1 + T1*B1;
      Nc = C1*Y1 - T1*A1;
  end;    
end;            

if revrsd,
   eval('[Nc,Dc] = reverse(Nc,Dc,''l'');','[Nc,Dc] = reverse(Nc,Dc,''r'');');
end;
Nc.v = var; Dc.v = var;
 
% References: [1] V. Kucera, Analysis and Design of Discrete Linear Control Systems,
%	      	  Academia, Prague, 1991
 
%end .. plqg
 