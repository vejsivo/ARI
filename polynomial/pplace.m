function [Y,X,D,E,degT] = pplace(B,A,C,arg3,arg4);
%PPLACE  Polynomial pole placement
%
% Given two polynomial matrices N and D and a vector of desired pole
% locations POLES the commands
%    [NC,DC] = PPLACE(N,D,POLES)
%    [NC,DC] = PPLACE(N,D,POLES,'l') 
% compute the polynomial matrix fraction representation C = -NC/DC 
% of a feedback controller for the system D\N that assigns the closed-loop 
% system poles to the location POLES. The multiplicity of the poles is
% increased if necessary. The resulting system may have real or complex 
% coefficients depending on whether or not the desired poles are self-conjugate.
%
% The commmand
%    [NC,DC] = PPLACE(N,D,POLES,'r') 
% computes the polynomial matrix fraction representation C = -DC\NC of such a 
% controller for the system P = N/D. 
%
% The commmand
%    [NC,DC,E,F,DEGT] = PPLACE(N,D,POLES), 
% possibly combined with the 'l' or 'r' options, additionally returns the 
% parametrization 
%	  K = -(NC+E*T)/(DC-F*T) 		('l' or empty option)
% or
%	  K = -(DC-T*F)\(NC+T*E)		('r' option)
% of all controllers that yield the same closed-loop poles. T is an 
% arbitrary polynomial matrix parameter of compatible size and of degree 
% bounded by degT. 
%
% NOTE: The user should be aware that for multi-input multi-output (MIMO) 
% systems just assigning pole locations need not be enough. In the MIMO case, 
% the desired behavior typically depends on the closed-loop invariant polynomials 
% rather than on the pole locations only. To assign the desired invariant 
% polynomials, put them into a diagonal matrix R of the same size as D and call 
%    [Nc,Dc] = pplace(N,D,R)
% or
%    [NC,DC,E,F,DEGT] = PPLACE(N,D,R).
% The 'l' or 'r' options may also be used in this case.  
%
% A tolerance TOL may be specified as an additional input argument.  

%	Author(s): M. Hromcik, M. Sebek 14-10-98
%	Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 1.0 $  $Date: 14-Oct-1998 10:28:34   $

global PGLOBAL;
eval('PGLOBAL.ZEROING;', 'painit;');
tol = PGLOBAL.ZEROING;

opt = 'l';
Cok = 0;

switch nargin
   
  case {0,1,2},
    error('Not enough input arguments.');
     
  case 3,
    if isa(C,'pol'), Cok = 1;
    end;
     
  case 4,
    if isa(C,'pol'), Cok = 1;
    end;
    if isnumeric(arg3) 
     tol = arg3;
    elseif ischar(arg3),
     opt = arg3;
    else 
     error('Invalid 4th argument.');
    end;
  
  case 5,
    if isa(C,'pol'), Cok = 1;
    end; 
    if isnumeric(arg3) 
     tol = arg3;
    elseif ischar(arg3),
     opt = arg3;
    else 
     error('Invalid 4th argument.');
    end;   
    if isnumeric(arg4) 
     tol = arg4;
    elseif ischar(arg4),
     opt = arg4;
    else 
     error('Invalid 5th argument.');
    end;
  
  otherwise
    error('Too many input arguments.');
end;	% switch  

if isempty(tol),
   tol = PGLOBAL.ZEROING;
elseif length(tol)~=1 | ~isreal(tol) | tol<0 | tol>1,
   error('Invalid tolerance.');
end;
if ~strcmp(opt,'r') & ~strcmp(opt,'l'),
   error('Invalid command option.');
end;

eval('A = pol(A); B = pol(B);', 'error(peel(lasterr));');

% A must be regular:
As = A.s; Bs = B.s;
rankA = rank(A, tol);  
if rankA < As(1),
  error('Denominator matrix is singular.');
end;

% Check dimensions
eval('[B,A] = testdnd(B,A,opt);', ...
   'error(peel(lasterr));');
   
% Check variable symbols and sampling periods
var = ''; h = 0;
eval('[var,h,B,A] = testvhnd(B,A);', ...
   'error(peel(lasterr));');
   
if ~Cok,
  % ..............
  % Compose C 1st:
  
  if ~isnumeric(C),
     error('Invalid 3rd argument.');
  end;
  poles = C(:).';  
  if isempty(poles),	% No poles prescribed
    error('No poles prescribed.');

  else
    switch opt,

      case 'l',
        [C,E,D, degT,bad_poles] = poles2c(A,B,var,h,poles,tol);
        eval('[X,Y] = axbyc(A,B,C,''miny'',tol);', 'error(peel(lasterr));');

      case 'r',
        [C,E,D, degT,bad_poles] = poles2c(A.',B.',var,h,poles,tol);  
        C = C.';
        E = E.'; D = D.'; degT = degT.';
        eval('[X,Y] = xaybc(A,B,C,''miny'',tol);', 'error(peel(lasterr));');

    end; 
    
    if any(any(isnan(X))) | any(any(isnan(Y))),
       error('Solution of polynomial equation failed;');
    end;
    if bad_poles, 
     warning(sprintf(['Some prescribed poles are not in the area of stability.\n', ...
     		'	 The closed loop system may be unstable.']));
    end; 		
   
  end;

else
  % .....................................
  % Solve diophantine equation with C(.):
  
  switch opt,

      case 'l',
        if nargout <=2,
          eval('[X,Y] = axbyc(A,B,C,''miny'',tol);', 'error(peel(lasterr));');
        else
          eval('[X,Y,E,D] = axbyc(A,B,C,''miny'',tol);', 'error(peel(lasterr));');
          E = -E;
          if nargout == 5,	% DEGT
            alpha = deg(A,'row');
            delta = deg(D,'col');
            xi = deg(C,'row') - alpha;
            xi = xi.';
            degT = zeros(As(1), size(D,1));
	         for i=1:size(D,1),
 	           degT(i,:) = xi - delta(i);
	         end;  
            degT(degT<0) = -Inf;
          end;
        end;

      case 'r',
        if nargout <= 2,
	       eval('[X,Y] = xaybc(A,B,C,''miny'',tol);', 'error(peel(lasterr));');
        else
          eval('[X,Y,E,D] = xaybc(A,B,C,''miny'',tol);', 'error(peel(lasterr));');
	       E = -E;
	       if nargout == 5,	% DEGT
            alpha = deg(A,'col');
            delta = deg(D,'row');
            xi = deg(C,'col') - alpha;
	         xi = xi.';
	         degT = zeros(As(1), size(D,1));
	         for i=1:size(D,1),
 	           degT(:,i) = xi - delta(i);
	         end;  
            degT(degT<0) = -Inf;
          end;  
	     end;
  end;
  
  if any(any(isnan(X))) | any(any(isnan(Y))),
     error('Solution of polynomial equation failed.');
  end;
  
  % Make X(z^-1) causal for D and Z^-1 variables:
  if any(strcmp(var, {'d', 'z^-1'})),
    degT = [];
    if rank(X{0}) < min(X.s),
      switch opt,
        case 'l',
          Z = A{0}\Y{0};
          if nargout <=2, 
            [E,D] = lmf2rmf(B,A,tol);
          end;
          X = X + E*Z;
          Y = Y - D*Z;  
        case 'r',
          Z = Y{0}/A{0};
          if nargout <=2, 
            [E,D] = rmf2lmf(B,A,tol);
          end;
          X = X + Z*E;
          Y = Y - Z*D;  
      end; 	%switch
    end;      
            
  else
    if rank(X) < size(X,1),
      warning(sprintf(['The denominator matrix of the controller is singular.\n',...
               '	 The desired controller probably does not exist.']));
    end; 
  end;            

end;

if ~isreal(X) | ~isreal(Y),
  warning('The compensator transfer function has complex coefficients.'); 
end;  

% ************************************************************* %
function [C,E,D,degT,bad_poles] = poles2c(A, B, var, h, poles, tol)
% POLES2C   Sub function for composing C of appropriate degrees
% from the prescribed roots in case of AX + BY = C - 'L' OPTION):

global PGLOBAL;

% Detect unstable poles:    
if strcmp(var,'s') | strcmp(var,'p'),
  bad_poles = any(real(poles) >= 0);
elseif strcmp(var,'z') | strcmp(var,'q')
  bad_poles = any(abs(poles) >= 1);
elseif strcmp(var,'z^-1') | strcmp(var,'d')
  bad_poles = any(abs(poles) <= 1);
else	% empty variable
  bad_poles = 0;  
end;  
    
As = A.s;

% A and B must be left prime
if ~isprime([A B], 'l'),
   error('Numerator and denominator are not coprime.');
  % G = gld(A,B,tol);
  % A = ldiv(A,G,tol);
  % B = ldiv(B,G,tol);  
end;

% A must be reduced:
if rank(lcoef(A,'row')) < As(1), 
   %error('A is not row reduced.'); 
   [A,pom,U] = rowred(A, tol);
   B = U*B;
end;
 
%% 
 
alpha = deg(A, 'row');
[E,D] = lmf2rmf(B,A);

delta = deg(D, 'col');
xi = delta(1) - 1;
	
gama = alpha + xi;

% Handle imaginary poles - add the conjugates:
poles_r = poles(~imag(poles));
poles_r = num2cell(poles_r);
poles_i = poles(~~imag(poles));

if ~isempty(poles_i),

  % Find conjugates:
  poles_i_p = sort(poles_i(imag(poles_i) > 0));
  poles_i_m = poles_i(imag(poles_i) < 0);
  [aux,ind] = sort(conj(poles_i_m));
  poles_i_m = poles_i_m(ind);	
  
  lenconj = min(length(poles_i_m), length(poles_i_p));
  poles_i = [ num2cell( [poles_i_p(1:lenconj); poles_i_m(1:lenconj)], 1 ), ...
              num2cell( poles_i_p(lenconj+1:end) ), ...
              num2cell( poles_i_m(lenconj+1:end) ) ];
              
  poles = [poles_r poles_i];
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % ADD conjugates:
  % poles_i = [poles_i; conj(poles_i)];
  % poles_i = num2cell(poles_i, 1);
  % poles = [poles_r poles_i];
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

else
  poles = poles_r;
end;  

rpt = ceil( sum(gama(gama>0))/length(poles) );
C = pol(zeros(As));
    
if rpt>1, 	% Not enough prescribed poles - repeat some of them:
   poles = repmat(poles, 1, rpt);
   if As(1)>0,
      for i=1:As(1)
         if gama(i) >= 0,	% else (i.e. if gama(i)<0), C(i,i) remains 0 
            len = 0;
            ptr = 1;
            polesi = {};
            while len<gama(i),
               polesi = {polesi{:}  poles{ptr}};
               len = len + length(poles{ptr});
               ptr = ptr+1;
            end; 
            polesi = cat(1, polesi{:});
            polesi = polesi.';
            Ci = root2pol(num2cell(polesi,2));
            Ci.v = var; Ci.h = h;
            if isempty(Ci), Ci = pol(1);
            end;
            C(i,i) = Ci;
            poles = poles(ptr:end);
         end;	%  if ~...  
      end;
   end;      
else	% Too many prescribed poles => make degCi > gama(i) for some i's        
   if As(1)>0,
      for i=1:As(1)
         len = 0;
         ptr = 1;
         polesi = {};
         while len<gama(i),
           polesi = {polesi{:}  poles{ptr}};
           len = len + length(poles{ptr});
           ptr = ptr+1;
         end; 
         polesi = cat(1, polesi{:});
         polesi = polesi.';
         Ci = root2pol(num2cell(polesi,2));
         Ci.v = var; Ci.h = h;
         if isempty(Ci), Ci = pol(1);
         end;
         C(i,i) = Ci;
         poles = poles(ptr:end);
     end;   
   end;      
   % Assign remaining roots:
   i = 0;
   while length(poles),
     i = mod(i, As(1))+1;
     rCi = root2pol(poles(1));
     rCi.v = var; rCi.h = h;
     C(i,i) = C(i,i) * rCi; 
     poles = poles(2:end);
   end;

end;  

% DEGT:
xi = deg(C,'row') - alpha;
xi = xi.';
degT = zeros(size(D,2), size(A,1));
D1 = size(D,1);
if D1>0,
   for i=1:D1,
     degT(i,:) = xi - delta(i);
   end;  
end;
degT(degT<0) = -Inf;
          
%end .. poles2c (sub function)          
          
% References: V. Kucera, P. Zagalak: Proper Solutions of Polynomial Equations, 
%	      submitted to ... , 1998.

%end .. pplace
