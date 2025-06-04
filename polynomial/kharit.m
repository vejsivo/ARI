function [stab,K1p,K2p,K3p,K4p,K1n,K2n,K3n,K4n] = kharit(Amin,Amax)
%KHARIT  Create the Kharitonov polynomials of an interval polynomial matrix
%
% The command
%    [STAB,K1P,K2P,K3P,K4P,K1N,K2N,K3N,K4N] = KHARIT(Amin,Amax) 
% creates eight fixed Kharitonov polynomials from the scalar interval 
% polynomial
%    A = [Amin,Amax]
% If  
%    A = [Amin{0},Amax{0}]+[Amin{1},Amax{1}]*v + ... +[Amin{d},Amax{d}]*v^d
% and
%    M1 = Amin{0}+Amin{1}*v+Amax{2}*v^2+Amax{3}*v^3+Amin{4}*v^4+Amin{5}*v^5+...
%    M2 = Amax{0}+Amax{1}*v+Amin{2}*v^2+Amin{3}*v^3+Amax{4}*v^4+Amax{5}*v^5+...
%    M3 = Amax{0}+Amin{1}*v+Amin{2}*v^2+Amax{3}*v^3+Amax{4}*v^4+Amin{5}*v^5+...
%    M4 = Amin{0}+Amax{1}*v+Amax{2}*v^2+Amin{3}*v^3+Amin{4}*v^4+Amax{5}*v^5+...
% then
%    K1P = real(M1) + j*imag(M4)         K1N = real(M4) + j*imag(M1)
%    K2P = real(M2) + j*imag(M3)         K2N = real(M3) + j*imag(M2)
%    K3P = real(M3) + j*imag(M1)         K3N = real(M2) + j*imag(M4)
%    K4P = real(M4) + j*imag(M2)         K3N = real(M1) + j*imag(M3)
%
% STAB is 1 if all eight Kharitonov polynomials are stable, and 0 otherwise.
%
% If Amin and Amax are matrices then the function works entry-wise.
%
% KHARIT serves to test the robust stability of continuous-time interval
% polynomials. It is useless for the discrete-time case.

%       Author(s):  S. Pejchova, M. Sebek 28-9-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 26-Oct-1998 11:13:34   $
%       $Revision: 3.0 $  $Date: 22-Aug-2000  S. Pejchova   $

global PGLOBAL;
eval('PGLOBAL.VARIABLE;', 'painit;');

ni=nargin;
narginchk(1,2);
% error(nargchk(1,2,ni));	%REMOVED IN NEW MATLABS
if ni==1,  Amax=Amin;  end;
vvk=PGLOBAL.VARIABLE;
eval('Amin=pol(Amin); Amax=pol(Amax);',...
   'error(peel(lasterr));');

[Asn1,Asn2,Adn]=size(Amin);
[Asx1,Asx2,Adx]=size(Amax);
Asn=[Asn1,Asn2]; Asx=[Asx1,Asx2];
Acn=Amin.c; Acx=Amax.c;  stab=0; cmpl=0;
if any(imag(Acn(:)))|any(imag(Acx(:))), cmpl=1; end
K1p=[]; K2p=[]; K3p=[]; K4p=[];K1n=[]; K2n=[]; K3n=[]; K4n=[];
vvn=Amin.var; vvx=Amax.var;
if isempty(vvn), vvn=vvx; end;
if isempty(vvx)|strcmp(vvn,vvx), 
   if ~isempty(vvn), vvk=vvn; end
else,
   warning('Inconsistent variables.');
end;

I1=[strmatch(vvk,{'z^-1';'d';'z';'q'},'exact')];
if ~isempty(I1),
   warning(['Discrete-time polynomial: robust stability can not ',...
            'be tested by Kharitonov polynomials.']);
end;
LAmin=lcoef(Amin,'ent'); LAmax=lcoef(Amax,'ent');   
if any(Asn - Asx),
   if all(Asn==1),
      Acn = repmat(Acn,Asx1,Asx2);
      Asn1=Asx1; Asn2=Asx2;
      LAmin=LAmin*ones(Asn1,Asn2);
   elseif all(Asx==1),
      Acx = repmat(Acx,Asn1,Asn2);
      LAmax=LAmax*ones(Asn1,Asn2);
   else,
      error('Matrices not of the same dimensions.');
   end;
end;
if isempty(Adn)|isempty(Adx), return;
elseif isinf(Adn) & isinf(Adx), 
   K1p=pol(zeros(Asn1,Asn2)); K2p=K1p; K3p=K1p; K4p=K1p;
   K1n=K1p; K2n=K1p; K3n=K1p; K4n=K1p; return;
elseif isinf(Adn), Acn=zeros(Asn1,Asn2); Adn=0;   
elseif isinf(Adx), Acx=zeros(Asn1,Asn2); Adx=0;
end;

if Adn~=Adx,
   dd=Adx-Adn; Adx=max(Adx,Adn);
   Acn=cat(3,Acn,zeros(Asn1,Asn2,dd));
   Acx=cat(3,Acx,zeros(Asn1,Asn2,-dd));
end;
ddn1=sign(real(LAmin)); ddx1=sign(real(LAmax));
ddn2=sign(imag(LAmin)); ddx2=sign(imag(LAmax));
dd1=(ddx1-ddn1)&((ddx2-ddn2)|(~ddn2));
dd2=(~ddn1)&(ddx2-ddn2);
if any(dd1(:))|any(dd2(:)),
   warning(['Degree drop: robust stability can not ',...
            'be tested by Kharitonov polynomials.']);
end;
C1=zeros(Asn1,Asn2,Adx+1);  C2=C1; C3=C1; C4=C1;
kxx=[1 1 0 0; 0 0 1 1; 0 1 1 0; 1 0 0 1];
for ii=1:4,
    C1(:,:,ii:4:end) = kxx(1,ii)*Acn(:,:,ii:4:end) + ...
                       kxx(2,ii)*Acx(:,:,ii:4:end);
    C2(:,:,ii:4:end) = kxx(2,ii)*Acn(:,:,ii:4:end) + ...
                       kxx(1,ii)*Acx(:,:,ii:4:end);
    C3(:,:,ii:4:end) = kxx(3,ii)*Acn(:,:,ii:4:end) + ...
                       kxx(4,ii)*Acx(:,:,ii:4:end);
    C4(:,:,ii:4:end) = kxx(4,ii)*Acn(:,:,ii:4:end) + ...
                       kxx(3,ii)*Acx(:,:,ii:4:end);
end;
K1p=real(C1)+j*imag(C4);  K1p=pol(K1p(:,:),Adx,vvk);
K2p=real(C2)+j*imag(C3);  K2p=pol(K2p(:,:),Adx,vvk);
K3p=real(C3)+j*imag(C1);  K3p=pol(K3p(:,:),Adx,vvk);
K4p=real(C4)+j*imag(C2);  K4p=pol(K4p(:,:),Adx,vvk);

K1n=real(C4)+j*imag(C1);  K1n=pol(K1n(:,:),Adx,vvk);
K2n=real(C3)+j*imag(C2);  K2n=pol(K2n(:,:),Adx,vvk);
K3n=real(C2)+j*imag(C4);  K3n=pol(K3n(:,:),Adx,vvk);
K4n=real(C1)+j*imag(C3);  K4n=pol(K4n(:,:),Adx,vvk);

stab=zeros(Asn1,Asn2); stab1=stab;
for ii=1:Asn1, for jj=1:Asn2,
      stab(ii,jj)=(isstable(K1p(ii,jj)))&(isstable(K2p(ii,jj)))&...
         (isstable(K3p(ii,jj)))&(isstable(K4p(ii,jj)));
      stab1(ii,jj)=(isstable(K1n(ii,jj)))&(isstable(K2n(ii,jj)))&...
         (isstable(K3n(ii,jj)))&(isstable(K4n(ii,jj)));     
   end;   end;
if cmpl, stab=stab&stab1; end

%end .. kharit
