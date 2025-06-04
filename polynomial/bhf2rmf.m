function [N,D] = bhf2rmf(F,G,H,bsizes,tol)
%BHF2RMF  Convert a block Hessenber form
%         into a right coprime polynomial matrix fraction
%
% The function
%    [N,D] = BHF2RMF(F,G,H,bsizes)
% obtains a right coprime matrix fraction description of a minimal 
% state space realization (A,B,C) in upper block Hessenberg form. 
%
% For bringing any observable state space realization into
% the appropriate upper block Hessenberg form,
% see BHF, FRAC/BHF.

% Author: R.C.W. Strijbos, November 13, 1998.
% Copyright 1998 by Polyx, Ltd.
% $ Revision 3.0 $  $ Date 09-Dec-2000  J.Jezek $
%                   $ Date 24-Feb-2003  J.Jezek $

% References: R.C.W. Strijbos, Calculation of Right Matrix Fraction
% Descriptions; an Algorithm, Proceedings of the 4th IEEE Mediterranean
% Symposium on New Directions in Control and Automation, Maleme,
% Krete, Greece, June 10-13, 1996, pp. 535-540.

% Initialisation of variables

global PGLOBAL;
eval('PGLOBAL.ZEROING;', 'painit;');

ni = nargin;
if ni<4,
   error('Not enough input arguments.');
end;
if ni==5 & ~isempty(tol),
   if ~isa(tol,'double') | length(tol)~=1 | ...
      ~isreal(tol) | tol<0 | tol>1,
      error('Invalid tolerance.');
   end;
else
   tol = PGLOBAL.ZEROING;
end;

[rG,cG]=size(G);
[rH,cH]=size(H);
mu=length(bsizes);
[Gtemp,T,rankG]=hhr(G);

% Begin Main Loop Here 

V=zeros(rG,cG*mu);
cumbsizes=cumsum(bsizes);
Vstart=bsizes(mu);
Vtemp = [eye(bsizes(mu)) zeros(bsizes(mu),cG-bsizes(mu))];	% Assign V_mu
iter=1;
while iter <= mu
   dmu = mu -iter;
   if dmu == 0
      rowstart = 1;
   else
      rowstart = cumbsizes(dmu)+1;
   end
   rowend = cumbsizes(dmu+1);
   V(rowstart:rowend,1:cG*iter) = Vtemp;
   RHS = [zeros(bsizes(dmu+1),cG) Vtemp];
   RHS = pol(RHS,iter);
   for k = 1 : iter
      if k<mu
	 colstart = cumbsizes(mu-k)+1;
      else
	 colstart = 1;
      end
      colend = cumbsizes(mu+1-k);
      Ftemp = F(rowstart:rowend,colstart:colend);
      Vtemp = pol(V(colstart:colend,1:cG*iter),iter-1);
      RHS = RHS - Ftemp * Vtemp;
   end
   RHS = RHS.c;
   if dmu ~= 0
      if dmu == 1
         colstart = 1;
      else
         colstart = cumbsizes(dmu-1)+1;
      end
      colend = cumbsizes(dmu);
      offset = bsizes(dmu)-bsizes(dmu+1);
      Ftemp = F(rowstart:rowend,colstart+offset:colend);
   else
      offset = cG-rankG;
      Ftemp = Gtemp(1:rankG,cG-rankG+1:cG);
   end
   Vtemp = [];
   if dmu ~= 0
      kstart = bsizes(dmu+1);
   else
      kstart = rankG;
   end
   for k = kstart : -1 : 1
      Vkk = zeros(1,cG*(1+iter));
      for kk = k+1 : kstart
	 Vkk = Vkk + Ftemp(k,kk) * Vtemp(kk-k,:);
      end
      Vtemp = [(RHS(k,:) - Vkk)/Ftemp(k,k) ; Vtemp];
   end
   for k = offset:-1:1
      Vkk = zeros(1,cG*(1+iter));
      Vkk(1,Vstart+k)=1;
      Vtemp = [ Vkk ; Vtemp];
   end
   iter = iter + 1;
   Vstart=Vstart+offset;
end

%End main Loop
   
D = pol(Vtemp,mu);
D = T * D;
V = pol(V,mu-1);
N = H * V;

var = D.v;        %  J.Jezek  24-Feb-2003     
if isempty(var), var = N.v;
end;
if strcmp(var,'z^-1') | strcmp(var,'d'),
   [N,D] = reverse(N,D,'r',tol);
end;              %  J.Jezek  24-Feb-2003

Nc = N.coef; Nme = norm(Nc(:,:),'inf');
Ntol = tol*Nme;
Dc = D.coef; Dme = norm(Dc(:,:),'inf');
Dtol = tol*Dme;

D = pzer(D,Dtol);
N = pzer(N,Ntol);

%end .. bhf2rmf
