function [F,G,H,bsizes] = bhf(A,B,C,tol)
%BHF  Converts the realization {A,B,C} into block Hessenberg form
%
% The function 
%    [F,G,H,BSIZES] = BHF(A,B,C,TOL)
% converts a state space realization (A,B,C) into upper block Hessenberg 
% form. The controllable part of the system {A,B,C} is returned as the 
% system {F,G,H}, where F is an upper block Hessenberg matrix and G has 
% nonzero elements in its first M rows only, where M is the rank of B. 
% The vector BSIZES contains the sizes of the different blocks of the 
% matrix F.
%
% A tolerance TOL may be specified as an additional input argument.

% Author: R.C.W. Strijbos, November 13, 1998.
% Copyright 1998 by Polyx, Ltd.

% References: R.C.W. Strijbos, Calculation of Right Matrix Fraction
% Descriptions; an Algorithm, Proceedings of the 4th IEEE Mediterranean
% Symposium on New Directions in Control and Automation, Maleme,
% Krete, Greece, June 10-13, 1996, pp. 535-540.

if nargin<3,
   error('Not enough input arguments.');
end;
if ~isa(A,'double') | ~isa(B,'double') | ~isa(C,'double') | ...
      ndims(A)>2 | ndims(B)>2 | ndims(C)>2,
   error('Invalid 1st, 2nd or 3rd argument; must be constant.');
end;
[rA,cA] = size(A); [rB,cB] = size(B); [rC,cC] = size(C);
if rA ~= cA | rB ~= rA | cC ~= rA
   error('Matrices of inconsistent dimensions.');
end;

bsizes = [];
if isempty(A) | isempty(B) | isempty(C),
   F = A; G = B; H = C;
   return;
end;

Tot = [B A];
if nargin == 3 | isempty(tol),
   normTot = norm(Tot);
   if normTot >= eps*1e2
      tol = (rA+rB)*normTot*eps*1e4;
   else
      tol = eps*1e2;
   end
else
   if ~isa(tol,'double') | length(tol)~=1 | ...
         ~isreal(tol) | tol<0 | tol>1,
      error('Invalid tolerance.');
   end;
end;

iter=1;
sumbsizes=rank(B,tol);
rowstart=1;
colstart=1;
colstop=cB;
T=eye(rA);
while (sumbsizes <= rA)
   Btemp=Tot(rowstart:rB,colstart:colstop);
   [rBt,cBt]=size(Btemp);
   Atemp=Tot(rowstart:rA,colstop+1:cA+cB);
   [rAt,cAt]=size(Atemp);
   [Btemp,Tloc,rankBtemp]=hhl(Btemp,tol);
   if rankBtemp == 0
      break
   end
   bsizes(iter)=rankBtemp;
   sumbsizes=sum(bsizes);
   iter=iter+1;
   T(rA-rAt+1:rA,:)=Tloc*T(rA-rAt+1:rA,:);
   Atemp=T*A*T';
   Btemp=T*B;
   Tot = [Btemp Atemp];
   rowstart=rowstart+rankBtemp;
   colstart=colstart+rankBtemp;
   colstop=colstop+rankBtemp;
end
F=T*A*T';
G=T*B;
H=C*T';
if rankBtemp == 0
   F=F(1:sumbsizes,1:sumbsizes);
   G=G(1:sumbsizes,:);
   H=H(:,1:sumbsizes);
end
iter=iter-1;
while iter>1
   rightstart=bsizes(iter-1)-bsizes(iter);
   if rightstart > 0
      break;
   else
      iter=iter-1;
   end
end
if iter > 1
   T=eye(sumbsizes);
   cumbsizes=cumsum(bsizes);
   Tloc=eye(bsizes(iter));
   while iter>1
      rowstart=cumbsizes(iter-1)+1;
      rowstop=cumbsizes(iter);
      if iter>2
         colstart=cumbsizes(iter-2)+1;
      else
         colstart=1;
      end
      colstop=rowstart-1;
      Ftemp=Tloc'*F(rowstart:rowstop,colstart:colstop);
      [Ftemp,Tloc,rankBtemp]=hhr(Ftemp,tol);
      T(colstart:colstop,colstart:colstop)=Tloc;
      iter=iter-1;
   end 
   F=T'*F*T;
   G=T'*G;
   H=H*T;
end

%end .. bhf

