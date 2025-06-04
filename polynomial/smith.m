function [S,U,V,UI,VI] = smith(A,varargin)
%SMITH  Smith diagonal form of a polynomial matrix
% 
% The commmands
%    [S,U,V,UI,VI] = SMITH(A,TOL)
%    [S,U,V,UI,VI] = SMITH(A,'elo',TOL)
% compute the Smith form of a polynomial matrix A by elementary operations.
%
% The commmand
%    [S,U,V] = SMITH(A,'moca',TOL)
% computes the Smith form of A by a 'Monte Carlo' method.
%
% The command
%    S = SMITH(A,'zeros',TOL)    
% returns the Smith diagonal form of a non-singular polynomial matrix A
% based on the computation of the zeros of A. The macro first computes
% the zeros of A with the macro ROOTS and then evaluates their algebraic 
% and geometric multiplicities with the macro CHARACT. The option 'zeros' 
% may be used ONLY FOR SQUARE MATRICES.
%
% If
%    A = A{0} + A{1}*v + A{2}*v^2 + ... + A{d}*v^d
% then its Smith diagonal form is
%                [ a11   0  ... ... ...  0  |
%                |  0   a22 ... ... ...  0  |
%    S = U*A*V = | ...  ... ... ... ... ... |
%                |  0    0  ... arr ... ... |
%                | ...  ... ... ...  0  ... |
%                |  0.  ... ... ... ...  0  ]
% The rank of A is r and the polynomial a(k,k) divides a(k+1,k+1)
% for k = 1, 2, .., r-1. U and V are unimodular transformation matrices 
% and their inverse matrices are UI = U^(-1) and VI = V^(-1).
%
% The 'Monte Carlo' method (option 'moca') does not return UI and VI.
% The method based on zeros (option 'zeros') does not return U, V, UI and VI.
%
% A tolerance TOL may be specified as an additional input argument.

%       Author(s):  S. Pejchova, D. Henrion 27-10-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 23-Sep-1999 15:05:34   $

global PGLOBAL;
eval('PGLOBAL.ZEROING;', 'painit;');
zr=PGLOBAL.ZEROING;
ni=nargin; method='elo'; lp=0; tol=[];
narginchk(1,4);
% error(nargchk(1,4,ni));	%REMOVED IN NEW MATLABS
eval('A=pol(A);','error(peel(lasterr));');
if ni>1,
   for ii=1:ni-1,
      argm=varargin{ii};
      if ischar(argm),
         if strcmp(argm,'moca')|strcmp(argm,'elo')|strcmp(argm,'zeros'),
            method=argm;
         elseif strcmp(argm,'second'), lp=10;
         else, error('Invalid command option.');
         end;
      elseif isa(argm,'double')& ~any(size(argm)-[1 1])& ...
            isreal(argm) & argm>=0 & argm<=1,
         tol=argm;
      else, error('Invalid tolerance.');
      end;
   end;
end;

[rS,cS,Sd]=size(A); ma=min(rS,cS); degofmat=Sd*ma;
U=pol(eye(rS));  UI=U;
V=pol(eye(cS));  VI=V;
if (isempty(Sd) | isinf(Sd)), S=A; return; end;
NormS=norm(A);

switch method,
case 'elo',
   S=A*(1/NormS);
   if isempty(tol),
      tol=(max(size(A{:})))*zr/10;
   else,
      tol=tol/NormS;
   end
   [DS,LS]=deg(S,'ent');
   loop=6*ma*(Sd+1);  lu=0;
   if ma >= 1
      ii=1;
      while ii<=ma                         % reduction of ii-row and ii-column
          xcol=0;, xrow=0;
          if max(DS(ii,ii+1:cS))>=0, xcol=1; end;
          if max(DS(ii+1:rS,ii))>=0, xrow=1; end;

          while (xcol | xrow) & (lu<loop)
             TU=eye(rS); TV=eye(cS);
             Min=min(min(abs(DS(ii:rS,ii:cS))));
             LSm=zeros(rS-ii+1,cS-ii+1);
             LSm(DS(ii:rS,ii:cS)==Min)=1;
             LSm=[zeros(ii-1,cS); zeros(rS-ii+1,ii-1),LSm];
             [Mi,rM]=max(abs(LS.*LSm)); [Mx,kj]=max(Mi);
             ki=rM(kj);                  % (ki,kj) - pivot
             if rS==1, kj=ki; ki=1; end;
             if (ki==ii) & (kj==ii) & ~isinf(DS(ii,ii)),
                TU=pol(TU); TUI=TU; TV=pol(TV); TVI=TV;
                ar=S(ii,ii);
                if rS > ii,
                   br=rdiv(S(ii+1:rS,ii),ar,tol);
                   TU(ii+1:rS,ii)=-br; TUI(ii+1:rS,ii)=br;
                end;
                if cS > ii,
                   br=rdiv(S(ii,ii+1:cS),ar,tol);
                   TV(ii,ii+1:cS)=-br; TVI(ii,ii+1:cS)=br;
                end;
                if any([any(isnan(TU)),any(isnan(TV))]),
                   error('The reduction failed. Try the option ''moca'' or ''zeros''.');
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
             S=mtimes(mtimes(TU,S,tol),TV,tol);
             U=mtimes(TU,U,tol);  V=mtimes(V,TV,tol);
             UI=mtimes(UI,TUI,tol);  VI=mtimes(TVI,VI,tol);
             [DS,LS]=deg(S,'ent');
             if (ki==ii) & (kj==ii) & ~isinf(DS(ii,ii)),
                Sd=deg(S); Sc=S{:};
                if cS>ii, Sc(ii,DS(ii,ii)*cS+ii+1:end)=0; end;
                if rS>ii, Sc(ii+1:rS,DS(ii,ii)*cS+ii:cS:end)=0; end;
                S=pol(Sc,Sd);
                [DS,LS]=deg(S,'ent');
             end  %(if k==ii)
             lu=lu+1; xrow=0; xcol=0;
             if max(DS(ii+1:rS,ii))>=0, xrow=1; end;
             if max(DS(ii,ii+1:cS))>=0, xcol=1; end;
          end  %(while lu<loop)
          if ii > 1
             u=S(ii-1,ii-1); w=S(ii,ii);
             ur=fliplr(u{:}); wr=fliplr(w{:});
             if norm(ur,inf) > tol
                [ddr,rrr]=deconv(wr,ur);
                if norm(rrr,inf) > tol,
                   g=[]; W=[];
                   eval('[g,W]=gld(u,w,''syl'',tol);',...
                      'error(lasterr);');
                   if isempty(g)|isempty(W),
                      error('The reduction failed. Try the option ''moca'' or ''zeros''.');
                   end;
                   p1=W(1,1); q1=W(2,1);
                   TU=pol(eye(rS)); TV=pol(eye(cS)); TUI=TU; TVI=TV;
                   TU([ii-1,ii],[ii-1,ii])=[p1 1; 1 0];
                   TUI([ii-1,ii],[ii-1,ii])=[0 1; 1 -p1];
                   TV(ii,ii-1)=q1;  TVI(ii,ii-1)=-q1;
                   S=mtimes(mtimes(TU,S,tol),TV,tol);
                   U=mtimes(TU,U,tol);  V=mtimes(V,TV,tol);
                   UI=mtimes(UI,TUI,tol);  VI=mtimes(TVI,VI,tol);             
                   S(ii-1,ii-1)=g;
                   [DS,LS]=deg(S,'ent');
                   ii=ii-2;
                end %(if norm(rrr))
             end %(if norm(ur)
          end %(if ii>1)
          ii=ii+1;
      end  %(while ii<=ma)
   end  %(if ma>1)
   if lu>=loop
      warning('Reduction interrupted.');  lp=1;
   end
   S=S*NormS; tol=tol*NormS;
   LS=lcoef(S,'ent'); dL=diag(LS);
   dL(find(dL==0))=1; dL=[dL; ones(cS-length(dL),1)];
   TV=diag(1./dL);
   S=mtimes(S,TV,tol);
   V=mtimes(V,TV,tol);
   VI=mtimes(diag(dL),VI,tol);
   
case 'moca',
   if isempty(tol), tol=zr; end;
   err_str=['Reduction failed. Try to modify tolerance or use option ''elo'' or ''zeros''.'];
   eval('[T,V]=hermite(A,tol); [S,U]=hermite(T.'',tol);', 'error(err_str);');
   U=U.';
   %%%%%% unsuccesful: trying third turn
   if norm((S.*eye(cS,rS))~=S)
      [S,W]=hermite(S,tol);
      U=mtimes((W.'),U,tol);
   end
   %%%%%% unsuccesful: trying again with transpose
   if norm((S.*eye(cS,rS))~=S)
      [S,W]=hermite(S,tol);
      U=mtimes((W.'),U,tol);
   end;
   S=S.'; 
   if norm((S.*eye(rS,cS))~=S),
      lp=lp+1;
   elseif ma>1,
      for ii=ma:-1:2,
         u=S(ii-1,ii-1); w=S(ii,ii);
         ur=fliplr(u{:}); wr=fliplr(w{:});
         if norm(ur,inf) > tol
            [ddr,rrr]=deconv(wr,ur);
            if norm(rrr,inf) > tol, lp=lp+1; break; end
         end; 
      end;
      if lp==10, lp=0;  end;
   end;
   if lp==1,
      W1=prand(0,'uni','int',rS);  W2=prand(0,'uni','int',cS);
      A=mtimes(mtimes(W1,A,tol),W2,tol);
      eval('[S,U,V]=smith(A,tol,''moca'',''second'');','error(lasterr);');
      U=mtimes(U,W1,tol); V=mtimes(W2,V,tol);
   elseif lp>1,
     error('No solution by the Monte Carlo method in this run. Try again or use the option ''elo'' or ''zeros''.');
   end;
   UI=[]; VI=[];
case 'zeros',
   if rS ~= cS,
      error('The input polynomial matrix must be square. Try the option ''moca'' or ''elo''.');
   end;
   if isempty(tol), tol=zr; end;
   U=[]; V=[]; UI=[]; VI=[];
   
   % Compute the zeros of A

   R = roots(A, tol);
   m = length(R);

   if m < 1,

      % Input polynomial matrix has no zeros
      S = pol(eye(rS));
   else

      % Sort the zeros by increasing real part
      [void, i] = sort(real(R)); R = R(i);

      % Extract distinct zeros
      m = length(R);
      R = real(R) + sqrt(-1) * sign(imag(R)) .* imag(R);
      Z = []; cpxZ = [];

      i = 1;
      while i <= m,

        z = R(i);
        cpx = (abs(imag(z)) > tol);
        j = sum(abs(R-z) < tol);

        if cpx & rem(j, 2), % complex zero with odd multiplicity
           error('Invalid complex zeros. Try the option ''moca'' or ''elo''.');
        end;

        Z = [Z z];
        cpxZ = [cpxZ cpx];

        % Next zero
        i = i + j;
      end;

      nz = length(Z); % number of distinct zeros

      % Build Smith form
      factors = pol(ones(1, rS));
      for i = 1:nz,
        z = Z(i);

        % Compute algebraic and geometric multiplicity
        % Increase tolerance until characteristic structure is found
        alg = 0; it = 0; newtol = tol;
        while (alg == 0) & (it < 10),
          [V, alg, geo, pow] = charact(A, z, newtol);
          it = it + 1; newtol = newtol * 10;
        end;

        if it == 10,
          error('Inaccurate zero computation. Try the option ''moca'' or ''elo''.');
        end;
        for j = 1:geo,
          k = rS-j+1;
          if (k < 1) | (k > rS),
            error('Invalid geometric multiplicity. Try the option ''moca'' or ''elo''.');
          end;

          if cpxZ(i),
            a = real(z); b = imag(z);
            factors(k) = factors(k) * pol([a^2+b^2 -2*a 1], 2) ^ pow(j);
          else
            factors(k) = factors(k) * pol([-z 1], 1) ^ pow(j);
          end; % complex
        end; % geom. mult.
      end; % dist. zeros
      S = diag(factors);
   end; % nb. zeros < 1
   
end;

if ~lp & ~isempty(U) & ~isempty(V),
  Res=(norm(minus(S,mtimes(mtimes(U,A,tol),V,tol),tol)))/NormS;
  if Res > tol*1e4,
     warning(sprintf('The relative residue of calculation is  %g',Res));
  elseif sum(sum(deg(S,'ent'))) > degofmat,
     warning('The resulting degree may not be correct!');
  end
end;

%end .. smith
