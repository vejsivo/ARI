function [D,U,V,UI,VI] = pdg(A,tol)
%PDG  Diagonalization of a polynomial matrix
% 
% The command
%    [D,U,V,UI,VI] = PDG(A,TOL)
% converts the polynomial matrix A by elementary operations
% to a diagonal matrix D so that U*A*V = D. U and V are 
% unimodular transformation matrices. The inverse matrices 
% UI=U^(-1) and VI=V^(-1) are also computed.
%
% A tolerance TOL may be specified as an additional input argument.

%       Author(s):  S. Pejchova 04-10-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 1.0 $  $Date: 20-Nov-1998 09:40:34   $
%       Modified by J. Jezek, 24-Aug-2001, variable, sampling per,

global PGLOBAL;
eval('PGLOBAL.ZEROING;', 'painit;');
zr=PGLOBAL.ZEROING;

ni=nargin;
narginchk(1,2);
% error(nargchk(1,2,ni));	%REMOVED IN NEW MATLABS
eval('A=pol(A);','error(peel(lasterr));');

var = A.v; per = A.h;
[rD,cD,Dd]=size(A); ma=min(rD,cD); degofmat=Dd*ma;
U=pol(eye(rD));  UI=U;
V=pol(eye(cD));  VI=V;
if (isempty(Dd) | isinf(Dd)), D=A; return; end;
NormD=norm(A);
D=A*(1/NormD);
if ni==2 & ~isempty(tol)
   if ~isa(tol,'double')| any(size(tol)-[1 1]) | ...
         ~isreal(tol) | tol<0 | tol>1,
       error('Invalid tolerance.');
   end;
   tol=tol/NormD;
else,
   tol=(max(size(A{:})))*zr/10;
end;

[DD,LD]=deg(D,'ent');
loop=6*ma*(Dd+1);  lu=0;
if ma >= 1
   ii=1;
   while ii<=ma                          % reduction of ii-row and ii-column
          xcol=0;, xrow=0;
          if max(DD(ii,ii+1:cD))>=0, xcol=1; end;
          if max(DD(ii+1:rD,ii))>=0, xrow=1; end;

          while (xcol | xrow) & (lu<loop)
             TU=eye(rD); TV=eye(cD);
             Min=min(min(abs(DD(ii:rD,ii:cD))));
             LDm=zeros(rD-ii+1,cD-ii+1);
             LDm(DD(ii:rD,ii:cD)==Min)=1;
             LDm=[zeros(ii-1,cD); zeros(rD-ii+1,ii-1),LDm];
             [Mi,rM]=max(abs(LD.*LDm)); [Mx,kj]=max(Mi);
             ki=rM(kj);                  % (ki,kj) - pivot
             if rD==1, kj=ki; ki=1; end;            
             if (ki==ii) & (kj==ii) & ~isinf(DD(ii,ii)),
                TU=pol(TU); TUI=TU; TV=pol(TV); TVI=TV;
                ar=D(ii,ii);
                if rD > ii,
                   br=rdiv(D(ii+1:rD,ii),ar,tol);
                   TU(ii+1:rD,ii)=-br; TUI(ii+1:rD,ii)=br;
                end;
                if cD > ii,
                   br=rdiv(D(ii,ii+1:cD),ar,tol);
                   TV(ii,ii+1:cD)=-br; TVI(ii,ii+1:cD)=br;
                end;             
                 if any([any(isnan(TU)),any(isnan(TV))]),
                   error('The reduction failed. Try to modify the tolerance.');
                end;
             elseif (ki~=ii)|(kj~=ii),
                TUI=TU; TVI=TV;
                if ki~=ii, 
                   TU([ii,ki],[ii,ki])=[0 -1; 1 0];
                   TUI([ii,ki],[ii,ki])=[0 1; -1 0];
                end
                if kj~=ii
                   TV([ii,kj],[ii,kj])=[0 -1; 1 0];
                   TVI([ii,kj],[ii,kj])=[0 1; -1 0];
                end
             end;
             D=mtimes(mtimes(TU,D,tol),TV,tol);
             U=mtimes(TU,U,tol);  V=mtimes(V,TV,tol);
             UI=mtimes(UI,TUI,tol);  VI=mtimes(TVI,VI,tol);
             [DD,LD]=deg(D,'ent');
             if (ki==ii) & (kj==ii) & ~isinf(DD(ii,ii)),
                Dd=deg(D); Dc=D{:};
                if cD>ii, Dc(ii,DD(ii,ii)*cD+ii+1:end)=0; end;
                if rD>ii, Dc(ii+1:rD,DD(ii,ii)*cD+ii:cD:end)=0; end;
                D=pol(Dc,Dd);
                D.v=var; D.h=per;
                [DD,LD]=deg(D,'ent');
             end  %(if k==ii)
             lu=lu+1; xrow=0; xcol=0;
             if max(DD(ii+1:rD,ii))>=0, xrow=1; end;
             if max(DD(ii,ii+1:cD))>=0, xcol=1; end;
          end  %(while lu<loop)
          ii=ii+1;
   end  %(while ii<=ma)
end  %(if ma>1)
if lu>=loop
   warning('Reduction interrupted.');
end
D=D*NormD;  tol=tol*NormD;
Res=(norm(minus(D,mtimes(mtimes(U,A,tol),V,tol),tol)))/NormD;
if Res > tol*1e4
   warning(sprintf('The relative residue of calculation is  %g',Res));
elseif sum(sum(deg(D,'ent')))  > degofmat
   warning('The resulting degree may not be correct!');
end;

%end .. pdg
