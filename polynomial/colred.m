function [D,rk,U,Ui] = colred(A, varargin)
%COLRED   Column-reduced form of a polynomial matrix
%
% The command
%    [D,RK,U,UI] = COLRED(A[,METHOD][,TOL]) 
% brings the polynomial matrix A into column reduced form D = A*U, 
% where U is a unimodular matrix and UI its inverse. The columns 
% of D are arranged according to decreasing degrees.
%
% If the input matrix A does not have full column rank then it is
% reduced to the form D = [D~ 0], with D~ column reduced.
%
% The number of nonzero columns in D is returned as the scalar RK.
%     
% The default method 'ref' ,is by repeated extraction of factors. A
% second method 'bas' is based on the calculation of a minimal basis.
%
% A tolerance TOL may be specified as an additional input argument.
%
% See also: ROWRED.

%    Authors: R.C.W. Strijbos,  D. Henrion, December 11, 1998.
%    Copyright 1998 by Polyx, Ltd.
%    $ Revision 3.0 $  $ Date 30-May-2000  J. Jezek $
%                      $ Date 08-Aug-2001  J. Jezek  sampl per  $
%    Modified by D. Henrion, December 19, 2002.
  
global PGLOBAL;
eval('PGLOBAL.VERBOSE;', 'painit;');

if nargin<1
   error('Not enough input arguments.');
end;
eval('D = pol(A);', 'error(peel(lasterr));');
typeD = D.var; perD = D.h;

if any(any(isnan(A))) | any(any(isnan(A))),
   error('Polynoial is not finite.');
end;

% verbose level
verbose = strcmp(PGLOBAL.VERBOSE, 'yes');

[rD, cD] = size(D); degD = D.degree;

% Default options.
method = 'ref';
tol = PGLOBAL.ZEROING;

% Parse optional input arguments.
lv = length(varargin);
if lv > 0,
   for i = 1:lv,
      arg = varargin{i};
      if isa(arg,'char') 
         switch arg
         case 'ref'
            method = 'ref';
         case 'bas'
            method = 'bas';
         otherwise
            error('Invalid command option.');
         end
      elseif (isa(arg,'double') & ~any(size(arg)-[1 1]))
         tol = arg;
         if ~isreal(tol) | tol<0 | tol>1,
            error('Invalid tolerance.');
         end
      else
         error(['Invalid ',nth(i+1),' argument.'])
      end
   end
end

if nargout == 4
   TUi = 1;
else
   TUi = 0;
end

% relative tolerance for zeroing
minD = min(abs(nonzeros(D.coef)));
if ~isempty(minD),
   tolzero = tol*minD;
else
   tolzero = 0;
end;

% REPEATED FRACTION EXTRACTION METHOD
if strcmp(method,'ref')
   U=pol(eye(cD),0); 
   if TUi, Ui=U; end

   if cD > 1
      cdeg=deg(D,'col');
      while (1 > 0) 
         [DegD,Dlc] = deg(D,'col');
         [E,IDPR]=cef(Dlc',tolzero*100);
         nullL=nullref(E,IDPR)';
         if isempty(nullL) 
            %cancel extraction procedure
            break; 
         end                                
         nullL=nullL/norm(nullL,1);
         maxDegD=max(DegD);
         Eps = 0;
         while (maxDegD>=0 & Eps < tolzero)
            dact=find(DegD==maxDegD);
            if ~isempty(dact)
	            Eps=max(max(abs(nullL(dact,:))));
               if Eps < tolzero
                  maxDegD=maxDegD-1;
               end
            else
	            maxDegD=maxDegD-1;
            end
         end
         if Eps ==0 , break, end
         [pc,qc]=find(abs(nullL(dact,:))==Eps);
         e=nullL(:,qc(1));
         k=dact(pc(1)); 
         Tk= pol(e /e(k),typeD); Tk.h = perD;
         for i=1:cD
            if (i~=k) & ~isinf(DegD(i)) 
               degs=DegD(k)-DegD(i);
               if degs >=0
                  sTk=shift(Tk(i),degs);
                  sTk.v = typeD; sTk.h = perD;
                  Tk(i)=sTk;
               else 
                  Tk(i)=0;
               end
            end
         end
         D(:,k) = D*Tk;
         U(:,k)=U*Tk;
         if TUi,
            T = pol(eye(cD)); T(:,k) = -Tk; T(k,k) = 1;
            Ui = T * Ui;
         end
         [Dga,Dla]=deg(D,'col');
         if Dga(k)==DegD(k)
            D{DegD(k)}(:,k)=zeros(rD,1);
         end;
      end   
   end
end;
% SYLVESTER MATRIX METHOD
if strcmp(method,'bas')
   if cD > 1,

      % The algorithm relies upon the computation
      % of the minimal polynomial basis of the right null space of
      % polynomial matrix [s^k*A -I] for a sufficiently high integer k. 
      cdeg = deg(D, 'col'); cdeg(cdeg < 0) = 0;
      k = sum(cdeg) - min(cdeg) + 1;
    
      % D = s^k*D
      Coef = D.Coef; 
      D = pol([zeros(rD, cD*k) Coef(:,:)], D.degree+k);
      D.v = typeD; D.h = perD;
      M = [D -eye(rD)];
      Z = null(M, tol);
   
      if isempty(Z),
       error('Reduction failed. Modify the tolerance.');
      end;

      % Z = [U; D]
      U = Z(1:cD, :);
      D = Z(cD+1:cD+rD, :);

      % D = D / s^k
      degD = D.degree;
      if degD < k,
         error('Incorrect kernel degree. Modify the tolerance.');
      end;
      Coef = D.coef; Coef = Coef(:,:,1+k:1+degD);
      D = pol(Coef(:,:), degD-k);
      D.v = typeD; D.h = perD;
   else
      % a matrix with only one column is necessarily column reduced

      U = 1;
      D = A;
      rk = norm(A) > tolzero;
   end
end;
D = pzer(D,tolzero);
U = pzer(U,tolzero);
if TUi
   if strcmp(method,'bas')
      Ui = pzer(pol(eye(cD)/U),tolzero);
   else
      Ui = pzer(Ui,tolzero);
   end
end
%   Sort Column degrees in descending order
cdeg=deg(D,'col');
[cdeg,I]   = sort(-cdeg);
D   = D(:,I);
U   = U(:,I);
if TUi, Ui = Ui(I,:); end
rk   = sum(isfinite(cdeg));

%Residue   check
if  verbose 
   residue = norm(D-A*U)/norm(A);
   disp(sprintf('Colred: Relative residue %g',residue));
end

%end .. colred



