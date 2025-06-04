function P = prand(varargin)
%PRAND Create a polynomial matrix with random real coefficients
%
% The command
%     P = PRAND(DEGP,I,J)  
% generates a "random"  I-by-J polynomial matrix P of degree DEGP 
% with normally distributed coefficients. If J is missing then a 
% square I-by-I matrix is created. If also I is missing then the
% sizes of DEGP are used.
%
% If DEGP is an integer matrix then its entries specify the degrees 
% of the entries of P. If DEGP is a column (row) vector, then its
% entries specify the row (column) degrees of P.
%
%     P = PRAND(DEGP,I,J,OPTIONS)  
% To generate the degrees of the entries, rows or columns of P
% randomly (within {0,DEGP}), include OPTIONS among the input arguments.
% OPTIONS is one of the strings 'ENT', 'ROW'  or 'COL'. 
%
%     P = PRAND(DEGP,I,'MON')  
% The OPTIONS = 'MON' generates P a square I-by-I monic polynomial matrix 
% of degree DEGP. 
%
%     P = PRAND(DEGP,I,'UNI')  
% The OPTIONS = 'UNI' makes P a square I-by-I unimodular matrix 
% (such that det(P) = 1).
%
%     P = PRAND(DEGP,I,'POS'[,ZEROS_VECTOR])  
% The input argument OPTIONS = 'POS'  generates again P a square I-by-I  
% polynomial matrix but now DEGP means the required number of zeros 
% including multiplicities. Some zeros can be fixed a priory by the vector
% ZEROS_VECTOR when complex conjugate complex parts are added if necessary.
%
%     P = PRAND(DEGP,I,'STA'[,ZEROS_VECTOR])  
% OPTIONS = 'STA' has the same functionality as 'POS' but P is
% 'stable' polynomial matrix (see the macro ISSTABLE for a more information).
% The a priory fixed zeros in ZEROS_VECTOR need not necessary be stable.
%
% Any of these syntaxes may be combined with the string 'INT' to produce
% "small integer" coefficients.

%       Author(s):  S. Pejchova, M. Sebek 31-08-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 13-Apr-2000 11:02:11   $

global PGLOBAL;
eval('PGLOBAL.VARIABLE;', 'painit;');
na = nargin;
if na > 6, error('Too many input arguments.');
elseif nargout > 1, error('Too many output arguments.');
end
i_tr = 0; typeA = 'full'; se_db = 0;
degA = max([1,round(abs(10*randn(1,1)))]);
rA = max([1,round(4*rand(1,1))]); 
cA = max([1,round(4*rand(1,1))]); 
vvx=PGLOBAL.VARIABLE;
if na>0,
 for kk=1:na
   argm=varargin{kk};
   if isa(argm,'double')&(ndims(argm)==2),
      argcol=argm(:); argcol((argcol<0)&(isinf(argcol)))=0;
      if any(~isfinite(argcol)),
         error('Nan and +Inf not allowed for degree and size.');
      elseif (se_db<2)&((norm(argcol-round(argcol))>PGLOBAL.ZEROING)|(any(argcol<0))),
         error('The degree and size must be nonnegative integers.');
      end;
      if ~se_db,
         degA=argm; se_db=1; rA=max([1 1;size(degA)]); cA=rA(2); rA=rA(1);
      elseif (isempty(argm))|(any(isinf(argm))),
            error('The size or the vector of zeros must be a nonempty and finite number.'); 
      elseif se_db==1,
         rA=argcol(1); cA=rA; se_db=2;
      elseif se_db==2,
         cA=argcol; se_db=3;
      end;
   elseif ischar(argm),
      switch argm,
       case {'ent','row','col','uni','sta','pos','mon'},
         typeA=argm;
       case {'s','d','p','q','z','z^-1','zi'},
          vvx=argm;
       case 'int', i_tr=1;
       otherwise,
         error('Invalid command option.');
      end;
   else,
      error(['Invalid ',nth(kk),' argument.']);
   end;
 end;
end;
[rd,cd]=size(degA);
if strcmp(typeA,'uni')|strcmp(typeA,'sta')|strcmp(typeA,'pos')|strcmp(typeA,'mon'),
   argcol=cA;, cA=rA;    
else,
   cA=round(abs(cA(1)));
   rd=min(rd,rA); cd=min(cd,cA); 
   degA=degA(1:rd,1:cd);
end;
maxdA=max(max(degA));
if (rA==0)|(cA==0)|(maxdA<0),
   P=pol(zeros(rA,cA));
   return;
elseif isempty(degA),
   error('Empty degree for nonempty polynomial matrix.');
elseif (rA==1)&(cA==1),
   switch typeA,
    case{'ent','row','col','full'},
      P=randn(rA,cA*(degA+1));
      if i_tr, P=round(5*P); end;
      P=pol(P,degA,vvx);  return;
    case 'uni', P=pol(1); return;
   end
end;
switch typeA,
 case 'uni',											% "unimodular" polynomial matrix
    no2=floor(maxdA/2);  no1=maxdA-no2;
    no4=1e3; no3=0.001; no5=0;
    if (rA*maxdA>10)&(rA*maxdA<=15), no4=1e2; no3=0.01;
    elseif rA*maxdA>15, no4=10; no3=0.1;
    end;
    if ((rA*maxdA>19) & ~i_tr) | (i_tr & (rA*maxdA>80)),
       no5=1; no1=maxdA; no2=no1-1; 
    end;
    A1=[zeros(rA,rA*(no2+1)),triu(no3*round(no4*randn(rA*(no1-no2))),1)];
    A2=zeros(rA,rA*(no2+1));
    for s1=1:no2+1,
       A12=no3*round(no4*randn(rA));
       A1(:,rA*(s1-1)+1:rA*s1)=triu(A12,1);           % upper(right) triang.
       A2(:,rA*(s1-1)+1:rA*s1)=tril(A12,-1);          % lower(left) triang.
    end;
    A1=A1+eye(rA,rA*(no1+1)); 
    A2=A2+eye(rA,rA*(no2+1));
    if i_tr, A1=round(A1); A2=round(A2); end;
    A2=pol(A2,no2,vvx);
    A1=pol(A1,no1,vvx);
    if no5,
       EP=eye(rA); EP1=EP(randperm(rA),:);EP2=EP(randperm(rA),:);
       EP2(1,:)=det(EP1)*det(EP2)*EP2(1,:);
       P=EP1*A1*EP2;
    else,
       P=A1*A2;
    end;
    return;
    
 case {'mon'},											% "monic" polynomial matrix
    P=randn(rA,rA*(maxdA));
    if i_tr, P=round(5*P); end;
    P=[P,eye(rA)];
    P=pol(P,maxdA,vvx);
    return;
    
 case {'sta','pos'},									% "stable" polynomial matrix
    maxdA=max([1,maxdA]); no1=0; nfre=0; nfco=0;
    fre=zeros(0,1); f_re=fre; f_im=fre;
    switch vvx,
    case {'z^-1','d','zi'}, no1=2;
    case {'q','z'}, no1=1;
    end;
    if se_db==3,
       fre=argcol(find(~imag(argcol))); nfre=length(fre);
       fco=find(imag(argcol));  nfco=length(fco);
       if nfco,
          f_co=sort(argcol(fco)); [f_x,ix]=sort(real(f_co));
          f_co=f_co(ix); f_re=real(f_co(1)); f_im=imag(f_co(1));
          ix=2;
          while ix<=nfco,
             if abs(f_co(ix)-conj(f_co(ix-1)))<eps, ix=ix+1; end;
             if ix < nfco+1,
                f_re=[f_re; real(f_co(ix))];
                f_im=[f_im; imag(f_co(ix))];
             end;
             ix=ix+1;
          end;
       end;
       nfco=length(f_re);
       maxdA=max([maxdA,nfre+2*nfco]);
    end;
    
    %Numbers of the different types of roots 
    %Probability of repeated roots is 0.05
    nrp=floor(sum(rand(maxdA,1)<0.05)/2);
    if se_db==3, nrp=0; end;
    
    %Probability of complex roots is 0.5 and imaginary part is bigger
    nco=floor(sum(rand(maxdA-2*nrp,1)<0.5)/2);
    if (se_db==3)&((nco<nfco) | ((maxdA-2*nco)<nfre)), nco=nfco; end;
    
    %Real roots
    nre=maxdA-2*nrp-2*nco;   
    
    rrp=-exp(randn(nrp,1));   rre=-2*exp(randn((nre-nfre),1));
    re=-exp(randn((nco-nfco),1));
    im=3*exp(randn((nco-nfco),1));
    if strcmp(typeA,'pos'),
       rrp=rrp.*sign(randn(nrp,1));
       rre=rre.*sign(randn((nre-nfre),1));
       re=re.*sign(randn((nco-nfco),1));
    end;
    
    rrx=[rrp;rrp;rre]; 

    if no1==1;
       rrx=(1+rrx)./(1-rrx);
       Zx=re+j*im; Zx=(1+Zx)./(1-Zx);
       re=real(Zx); im=imag(Zx);
    elseif no1==2,
       rrx=(1-rrx)./(1+rrx);
       Zx=re+j*im; Zx=(1-Zx)./(1+Zx);
       re=real(Zx); im=imag(Zx);
    end;
    if i_tr,
       if strcmp(typeA,'pos')&(no1==0 |no1==2),
          rrx=round(rrx); re=round(re); im=round(im); 
       else,
          switch no1,
          case 0,
             rrx=floor(rrx);  re=floor(re); im=round(im); 
          case 1,
             rrx=(fix(10*rrx))/10;  re=(fix(10*re))/10; im=(fix(10*im))/10; 
          case 2,
             rrx=round(rrx+(0.5*sign(rrx)));  
             re=round(re+(0.5*sign(re))); 
             im=round(im+(0.5*sign(im))); 
          end;
       end;
    end;
    if se_db==3,
       rrx=[rrx; fre];
       re=[re; f_re]; im=[im; f_im];
    end;

    nrx=length(rrx);
    rP=rA;
    if (rA==1) & (nrx+nco), rP=nrx+2*nco; end;
    P=pol(eye(rP),vvx);
    while nco | nrx,
       if nco,
          P1=eye(rP); Paux=zeros(rP);
          nco1=min([nco,(rP-mod(rP,2))/2]);
          nco=nco-nco1; re1=re(1:nco1); im1=im(1:nco1);
          d1=ones(nco1,1);
          if i_tr, [n1,d1]=rat(re1+i*im1); re1=real(n1); im1=imag(n1); end;
          Pr=[-re1';re1'];
          Pi=[im1'; zeros(1,nco1)]; Pi=Pi(:);
          Pi=diag(Pr(:))+diag(Pi(1:end-1),1)+diag(Pi(1:end-1),-1);
          Pr=[d1'; -d1']; Pr=diag(Pr(:));
          P1(1:2*nco1,1:2*nco1)=Pi;
          Paux(1:2*nco1,1:2*nco1)=Pr;
          P=P*pol([P1,Paux],1,vvx);
          re(1:nco1)=[]; im(1:nco1)=[];
       end;
       if nrx,
          P1=eye(rP); Paux=zeros(rP);
          nrx1=min([rP,nrx]);
          nrx=nrx-nrx1; n1=rrx(1:nrx1); d1=ones(nrx1,1);
          if i_tr, [n1,d1]=rat(rrx(1:nrx1)); end;
          idx=[rP-nrx1+1:rP];
          P1(idx,idx)=diag(-n1); Paux(idx,idx)=diag(d1);
          P=P*pol([P1,Paux],1,vvx);
          rrx(1:nrx1)=[];
       end;
    end;
    if rA==1, 
       P=det(P);
       if P{deg(P)}<0, P=-1*P; end;
    else,
        degPU=(rand(1,1)<0.5)+2*(rand(1,1)<0.1)+3*(rand(1,1)<0.01);
        if (maxdA+rA)>10, degPU=0; end;
        UP=prand(degPU,'uni','int',rA,vvx);
        P=P*(UP.');
     end;
     
    return;

 case 'full',
   if rd==1,
      DegP=[repmat(degA,rA,1),repmat(maxdA,rA,cA-cd)];
   elseif cd==1,
      DegP=[repmat(degA,1,cA); repmat(maxdA,rA-rd,cA)];
   else,
      DegP=repmat(maxdA,rA,cA);
      DegP(1:rd,1:cd)=degA;
   end;
 case 'row',										% row degrees "randomly"
   DegP=round((maxdA+1)*rand(rA,1))-1;
   DegP=repmat(DegP,1,cA);
   if (rd==1) & (cd~=1),
      DegrP=[repmat(degA,rA,1), repmat(maxdA,rA,cA-cd)];
      DegP=min(DegP,DegrP);
   end;
 case 'col',										% column degrees "randomly"
   DegP=round((maxdA+1)*rand(1,cA))-1;
   DegP=repmat(DegP,rA,1);
   if (cd==1) & (rd~=1),
      DegcP=[repmat(degA,1,cA); repmat(maxdA,rA-rd,cA)];
      DegP=min(DegP,DegrP);
   end;
 case 'ent',										% entries degrees "randomly"
   DegP=round((maxdA+1)*rand(rA,cA)); DegP=DegP-1;
end; %switch typeA

P=repmat(zeros(rA,cA),[1 1 (maxdA+1)]);
for ii=1:rA, for jj=1:cA
      if DegP(ii,jj)>=0,
         P(ii,jj,1:DegP(ii,jj)+1)=randn(1,1,(DegP(ii,jj)+1));
      end;
end; end;
if i_tr, P=round(5*P); end;
if maxdA >0, P=pol(P(:,:),maxdA,vvx);
else, P=pol(P(:,:));
end;

%end .. prand
