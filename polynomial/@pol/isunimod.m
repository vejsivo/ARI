function r = isunimod(P,arg2)
%ISUNIMOD Test if polynomial matrix is unimodular
%
% The command
%    R = ISUNIMOD(P[,TOL])
% returns 1 if the polynomial matrix P is unimodular and 0 if it is not.
%
% An optional tolerance TOL may be included. Its default value is the 
% global zeroing tolerance.
%
% See also POL/ISPRIME.

%	Author(s): M. Hromcik, H. Kwakernaak, M. Sebek 30-10-98
%	Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 19-Nov-1998 10:28:34   $
%                         $Date: 28-Mar-2003 warning    $

global PGLOBAL;

eval('tol = PGLOBAL.ZEROING;', 'painit; tol = PGLOBAL.ZEROING;'); 

switch nargin,
    case 0, error('Not enough input arguments.');
    case 1,
    case 2,
       if ~isempty(arg2),
          if isnumeric(arg2) & length(arg2)==1 & isreal(arg2) & ...
                arg2>=0 & arg2 <=1,
  	          tol = arg2;
          else
    	       error('Invalid tolerance.');
          end;
       end;    
end;		%switch

% Initializations
P = pol(P);
[m,n] = size(P);
if n ~= m
   error('Matrix is not square.')
end
   
% Test unimodularity:
detP = det(P,0);
detPc = detP.c;
  
if ne(detP, 0, tol) & all(isfinite(detPc(:))),		
  
    % detP is nonzero with tolerance TOL
    % ..................................
    r = all( abs(detPc(2:end)) < tol*abs(detPc(1)) );
    return;
else

    % First try the rank using numerically reliable method:
    switch nargin,
      case 1,
      	rankP = rank(P);
      case 2,
	rankP = rank(P,tol);
    end;		%switch

    if rankP < min(m,n),	% Singular matrix
      r = 0;
      warning('Matrix is singular to working precision.');
      return;
    end;  
  
    % P is badly scaled - is non singular and has zero determinant
    % Try to scale the input P to make its determinant 
    % "reasonable" - with unity leading or closing coefficient:
    % ............................................................    
    
    Plc = lcoef(P,'col');
    saved_w = warning; warning off;
    Plci = inv(Plc); warning(saved_w);
    
    if all(all(isfinite(Plci))),	% col reduced
    %disp col_red
      % .............................................
      % P is col reduced - l.c. of det(P) is det(Plc)
      % .............................................
      		
      P = P * Plci;	% det(P*Plci) is monic.
      detP = det(P,0);
      detPc = detP.c;
      r = all( abs(detPc(2:end)) < tol*abs(detPc(1)) );
      return;
      
    else
      Plr = lcoef(P,'row');
      saved_w = warning; warning off;
      Plri = inv(Plr); warning(saved_w);
      
      if all(all(isfinite(Plri))),		% row reduced
      %disp row_red
        % .............................................
        % P is row reduced - l.c. of det(P) is det(Prc)		
        % .............................................
        
        P = P * Plri;	% det(P*Plri) is monic
        detP = det(P,0);
        detPc = detP.c;
        r = all( abs(detPc(2:end)) < tol*abs(detPc(1)) );
        return;
     
      else 
        Pflip = P;
        Pflip.v = 'd';
        Pflip = Pflip';	
        Pcc = lcoef(Pflip,'col').';
        saved_w = warning; warning off;
        Pcci = inv(Pcc); warning(saved_w);
        
        if all(all(isfinite(Pcci))),	% back col reduced
        %disp back_col
          % ..................................................
          % P is back-col reduced - c.c. of det(P) is det(Pcc)
          % ..................................................
      		
          P = P * Pcci;		% det(P*Pcci) has unity closing coefficient
          detP = det(P,0);
          detPc = detP.c;
          r = all( abs(detPc(2:end)) < tol*abs(detPc(1)) );
          return;
     
        else
          Pcr = lcoef(Pflip,'row').';
          saved_w = warning; warning off;
          Pcri = inv(Pcr); warning(saved_w);
          
          if all(all(isfinite(Pcri))),		% back row reduced
          %disp back_row
            % ..................................................
            % P is back-row reduced - c.c. of det(P) is det(Pcr)		
            % ..................................................
          
            P = P * Pcri;	% det(P*Pcri) has unity closing coefficient
            detP = det(P,0);
            detPc = detP.c;
            r = all( abs(detPc(2:end)) < tol*abs(detPc(1)) );
            return;
  
          else
            % ..........................................................
            % P not reduced - do not try reduction (only little chance), 
            % return 0 and cause warning
            % ..........................................................
         
            % eval('[R,U] = colred(P);', 'error(''The matrix is badly scaled.'');');  
            % Rl = lcoef(R, 'col');
            % Rli = inv(Rl);
            % P = R * Rli;
            % detP = det(P,0);
            % detPc = detP.c;
            % r = all( abs(detPc(2:end)) < tol*abs(detPc(1)) );
            % return;
        
            warning('Matrix is badly scaled.');
            r = 0;
            
          end;					% ..end /back row reduced
        end;				% ..end  /back col reduced
      end;			% ..end /row reduced
    end;		% ..end /col reduced  
     
end;  		
 
%end .. @pol/isunimod
 
