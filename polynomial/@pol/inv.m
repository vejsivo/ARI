function M = inv(P,arg2,arg3)
%INV    Inverse of polynomial
%
% The command  M = INV(P)  returns the inverse of square polynomial
% matrix P. The result is a scalar-denominator fraction.
%
% An optional input argument TOL may specify a zeroing tolerance
% to be used instead of PGLOBAL.ZEROING.
%
% An optional input argument MET may specify the method for computing
% the adjoint. Methods are  'def' (default), 'int' .
%
% See also POL/ADJ, SDF/SDF.

%       Author:  J. Jezek  26-Jan-2000
%       Copyright(c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 21-Apr-2000 $
%                     $ Date 07-Nov-2000 $
%       $ Revision $  $ Date 07-Dec-2001 $ M.Hromcik - inv_local added at the end of file
%                     $ Date 30-Sep-2002 $
%                     $ Date 28-Feb-2003 $

global PGLOBAL;

tol = PGLOBAL.ZEROING; met = 'def';

if nargin==3 & ~isempty(arg3),
   if isa(arg3,'char'), met = arg3;
   elseif isa(arg3,'double'), tol = arg3;
   else error('Invalid 3rd argument.');
   end;
end;
if nargin>=2 & ~isempty(arg2),
   if isa(arg2,'char'), met = arg2;
   elseif isa(arg2,'double'), tol = arg2;
   else error('Invalid 2nd argument.');
   end;
end;

if length(tol)~=1 | ~isreal(tol) | tol<0 | tol>1,
   error('Invalid tolerance.');
end;

Adj = 0; Det = 0; M = 0;
eval('[Adj,Det] = inv_local(P,tol,met);','error(peel(lasterr));');
if Det==0,
   error('Determinant of matrix is zero.');
end;
M = sdf(Adj,Det);

if strcmp(PGLOBAL.COPRIME,'cop'),
   M = coprime(M,tol);
end;
if strcmp(PGLOBAL.REDUCE,'red'),
   M = reduce(M);
else
   M = smreduce(M);
end;
if strcmp(PGLOBAL.DEFRACT,'defr'),
   M = defract(M);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ad, dA] = inv_local(A, tol, met)
% INV  Inverse of a polynomial matrix as implemented in V 2.5

%	Author(s): M. Hromcik, M. Sebek 16-9-98
%	Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 1.0 $  $Date: 16-Sep-1998 10:28:34   $
%       Modified by S. Pejchova, February 19, 2001.
%                by J. Jezek, 24-Feb-2003, arg check

global PGLOBAL;

As = A.s; Ac = A.c;
rankA = rank(A,tol);

if rankA < As(1),	% singularity
  % ...............
  % singular matrix
  % ...............
  
  if As(1)-As(2),
    error('The matrix must be square.');
  end;
  error('Matrix is singular to working precision.');

elseif all(As==0),	
  % ...................
  % empty 0-by-0 matrix
  % ...................

  Ad = pol([]);
  dA = pol(1);

else			
    % ...................
    % non singular matrix
    % ...................
    Ad = 0; dA = 0;
    % ........................
    % Try scaling if A reduced
    % ........................
  
    Alc = lcoef(A,'col');
    warning off; Alci = inv(Alc); warning on;
    
    
    if rcond(Alci) > 1e-14,
    %if all(all(isfinite(Alci))),	% col reduced
      %disp col_red
      % .............................................
      % A is col reduced - l.c. of det(A) is det(Alc)
      % .............................................
      
      A = A * Alci;	
      eval('[Ad,dA] = adj(A,met,tol);','error(peel(lasterr));');
      Ad = Alci * Ad;
    else
      Alr = lcoef(A,'row');
      warning off; Alri = inv(Alr); warning on;

      if rcond(Alri) > 1e-14,
      %if all(all(isfinite(Alri))),		% row reduced
        %disp row_red
        % .............................................
        % A is row reduced - l.c. of det(A) is det(Arc)		
        % .............................................
        
        A = A * Alri;
        eval('[Ad,dA] = adj(A,met,tol);','error(peel(lasterr));');
        Ad = Alri * Ad;
      else 
        %error('cau');	
        Aflip = A;
        Aflip.v = 'd';
        Aflip = Aflip';
        Acc = lcoef(Aflip,'col');
        warning off; Acci = inv(Acc); warning on;				

        if rcond(Acci) > 1e-14,
        %if all(all(isfinite(Acci))),	% back col reduced
          %disp back_col
          % ..................................................
          % A is back-col reduced - c.c. of det(A) is det(Acc)
          % ..................................................
      		
          A = A * Acci;
          eval('[Ad,dA] = adj(A,met,tol);','error(peel(lasterr));');
          dc = dA.c;
          dA = pzer(dA, tol*max(abs(dc(:))));
          dc = dA.c;
          Ad = Acci * Ad;
          k = dc(end);
          dA.c = dc/k;
          Ad.c = Ad.c/k;
        else
          Acr = lcoef(Aflip,'row');
          warning off; Acri = inv(Acr); warning on;				
          
          if rcond(Acri) > 1e-14,
          %if all(all(isfinite(Acri))),		% back row reduced
          %disp back_row
            % ..................................................
            % A is back-row reduced - c.c. of det(A) is det(Acr)		
            % ..................................................
          
            A = A * Acri;
            eval('[Ad,dA] = adj(A,met,tol);','error(peel(lasterr));');
            dc = dA.c;
            dA = pzer(dA, tol*max(abs(dc(:))));
            dc = dA.c;
            Ad = Acri * Ad;  
            k = dc(end);
            dA.c = dc/k;
            Ad = Ad/k;
          else
            % ..........................................................
            % A not reduced - do not try reduction (only little chance), 
            % return NaN and cause warning
            % ..........................................................
         
            % eval('[R,U] = colred(A);', 'error(''The matrix is badly scaled.'');');  
            % Rl = lcoef(R, 'col');
            % Rli = inv(Rl);
            % Pom = R * Rli;
            % [APom, dA] = adj(Pom,met,tol);
            % Ad = U * Rli * APom;
            % Adc = Ad.c;
            % Ad = pzer(Ad, tol*max(abs(Adc(:))));
            % dc = dA.c;
            % k = dc(end);
            % dA.c = dA.c/k;
            % Ad.c = Ad.c/k;  
            
            eval('[Ad,dA] = adj(A,met,tol);','error(peel(lasterr));');            

            if dA==0, 
              	warning('Matrix is badly scaled.');
              	Ad = pol(NaN*eye(As));
              	dA = 1;
              	return;
            else				
    		% ....................................
    		% det(A) is nonzero - norm denominator
    		% ....................................
    
    		dc = dA.c;
    		dA = pzer(dA, tol*max(abs(dc(:))));
    		dc = dA.c;
    		k = dc(end);
    		dA.c = dc/k;
    		Ad.c = Ad.c/k;
   	    end; 	
         
          end;					% ..end /back row reduced
        end;				% ..end  /back col reduced
      end;			% ..end /row reduced
    end;		% ..end /col reduced 
end;			% ..end /singularity  
  
%end .. inv_local
  
%end .. @pol/inv
