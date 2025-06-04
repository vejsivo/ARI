function [D,L] = deg(A,string)
%DEG  Various degrees matrices of a polynomial matrix:
%
% D = DEG(A)       default, the same as DEG(A,'mat').
% D = DEG(A,'mat') returns the degree of polynomial matrix A.
% D = DEG(A,'ent') returns the matrix of degrees of the A entries.
% D = DEG(A,'row') returns the column vector of row degrees of A.
% D = DEG(A,'col') returns the row vector of column degrees of A.
% D = DEG(A,'dia') returns for a para-Hermitian polynomial matrix
% A(s)=A'(-s) a vector of half diagonal degrees of A.
%
% The second output argument L in all cases returns corresponding
% matrix or vector of coefficients, the same as first output argument
% in function LCOEF.
%
% See also LCOEF.

%       Author(s):  S. Pejchova, M. Sebek 13-3-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 22-Apr-1999 12:18:34   $

% Effect on other properties:
% D and L are standard Matlab matrices.

if nargin==1,
   string='mat';
elseif ~ischar(string),   
  error('Illegal input argument for method.');
end;
no=nargout;

Ac=A.c; Ad=A.d;  As1=A.s(1); As2=A.s(2);
L=zeros(As1,As2);

if isempty(Ad), D=[]; return; end;
if strcmp(string,'ent'),
   if As1==1, string='col';
   elseif As2==1, string='row'; end;
end;
switch string
 case 'mat',
      D=Ad;
      if (no==2)&(D>=0), L=Ac(:,:,D+1); end;

 case 'ent',
      if isinf(Ad),
         D=repmat(-Inf,[As1,As2]);
      else,
         Do=abs(Ac)&abs(Ac);
         if Ad>0,
            a = repmat(reshape(1:Ad+1,[1 1 Ad+1]),[As1 As2]);
            Do = Do.*a; Do = max(Do,[],3); 
         end;
         D = Do-1;          
         D(D<0)=-Inf;
         if no==2,
            L = Ac;
            if Ad > 0,
               Do=((repmat(Do,[1,1,Ad+1]))==a);
               L=sum((Do.*Ac),3);
            end;
         end;         
      end;

 case 'row',
      if isinf(Ad),
         D=repmat(-Inf,[As1,1]);
      else,
         Do=abs(Ac)&abs(Ac);
         if Ad>0,
            pD=permute(sum(Do,2),[1,3,2]);
            pD1=pD&pD; a=diag(1:Ad+1);
            Do=pD1*a;  
         end;
         Do = max(Do,[],2);
         D = Do-1;
         D(D<0)=-Inf; Do(Do==0)=1;
         if no==2,
            if Ad > 0,
              for ii=1:As1,  L(ii,:)=Ac(ii,:,Do(ii)); end;
            else,
               L = Ac;
            end;
         end;         
      end;
 case 'col',
      if isinf(Ad),
         D=repmat(-Inf,[1,As2]);
      else,
         Do=abs(Ac)&abs(Ac);
         if Ad>0, 
            pD=permute(sum(Do,1),[3,2,1]);
            pD1=pD&pD; a=diag(1:Ad+1);
            Do=a*pD1;  
         end;
         Do=max(Do,[],1);
         D = Do-1;
         D(D<0)=-Inf; Do(Do==0)=1;           
         if no==2,
            if Ad > 0,
               for ii=1:As2,  L(:,ii)=Ac(:,ii,Do(ii)); end;
            else,
               L = Ac;
            end;
         end;         
      end;
 case 'dia',
      global PGLOBAL;
      if As1~=As2,
         error('The input matrix is not square.');
      elseif norm(A-A','blk',inf) > PGLOBAL.ZEROING,
         error('The input matrix is not para-Hermitian');
      end;
      if isinf(Ad),
         D=repmat(-Inf,[As1,1]);
      else,
         Do=Ac&Ac;
         if Ad>0,
           a = repmat(reshape(1:Ad+1,[1 1 Ad+1]),[As1 As2]);      
           Do = Do.*a; Do = max(Do,[],3);
        end;
         D=(Do-1).*((Do-1)>=0);
         halfdiag = 1;
         for i = 1:As1, for j = i+1:As1,
             if 2*D(i,j) > D(i,i)+D(j,j), halfdiag = 0; end;
         end; end;
	
         if halfdiag,
            % if halfdiag = 1 diagonal degrees are half the
            % degrees of diagonal entries

            D = diag(D)/2;

         else
            warning('The half diagonal degrees are not unique');
            % computation of the smallest diagonal degrees sdiag(i)
            % such that D(i,j) <= di+dj
            ddiag = [D(1,1)/2 zeros(1,As1-1)];
            for i = 2:As1,
                ddiag(i) = max([D(i,1:i-1)-ddiag(1:i-1) D(i,i)/2]);
            end;
            D=ddiag';
         end; % if halfdiag
         if no==2,
            for i = 1:As1, for j = 1:As1,
                if D(i)+D(j) <= Ad,
                   b = Ac(i,j,(D(i)+D(j)+1));
                   L(i,j) = (-1)^D(i)*b;
                end;
            end; end;
         end;

      end;

 otherwise,
   error('Illegal input argument for method.');

end; %switch string

% end ../@pol/deg.m

