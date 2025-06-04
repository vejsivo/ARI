function [X,Y] = axby0(A,B,tol);
%AXBY0    Solver of polynomial equation  A*X + B*Y = 0
%
% For polynomial matrices A and B,
% number of rows of A == number of rows of B,
% matrix A square and nonsingular,
% the command
%     [X,Y] = AXBY0(A,B)
% returns polynomial matrices X and Y satisfying the equation
%     A*X + B*Y = 0 ,
% number of columns of X == number of columns of Y,
% matrix Y square and nonsingular.
% 
% A tolerance may be specified as an additional input argument.
% Its default value is the global zeroing tolerance.
%
% See also XAYB0.

%      Author:  J. Jezek  19-Nov-1999
%      Copyright(c) 1999 by Polyx Ltd.
%      $ Revision $  $ Date 24-May-2000 $
%                    $ Date 09-Dec-2000 $
%                    $ Date 22-Jul-2001 $

global PGLOBAL;
eval('PGLOBAL.ZEROING;', 'painit;');

ni = nargin;
if ni<2,
   error('Not enough input arguments.');
elseif ni==2 | isempty(tol),
   tol = PGLOBAL.ZEROING;
else
   if ~isa(tol,'double') | length(tol)~=1 | ...
         ~isreal(tol) | tol<0 | tol>1,
      error('Invalid tolerance.');
   end;
end;

eval('A = pol(A); B = pol(B);', ...
   'error(peel(lasterr));');
eval('[B,A] = testdnd(B,A,''l'');', ...
   'error(peel(lasterr));');

[tv,Xv,A,B] = testvp(A,B);
if tv==0,
   warning('Inconsistent variables.');
end;

[th,Xh,A,B] = testhp(A,B,Xv);
if th==0,
   warning('Inconsistent sampling periods.');
end;

[rA,cA] = size(A); [rB,cB] = size(B);
Ac = A.coef; Bc = B.coef;
Ame = norm(Ac(:,:),'inf'); Bme = norm(Bc(:,:),'inf');
if Ame~=0 & Bme~=0, kappa = Ame/Bme; B = B*kappa;
else kappa = 1;
end;
tolzero = tol*Ame;
XY = null(horzcat(A,B),tolzero);
cXY = size(XY,2);
if cXY<cB,
   error('axby0 system error: solution not found.');
elseif cXY>cB,
   XY = XY(:,1:cB);
end;
%X = pzer(XY(1:rB,:),tolzero);         Troubles, the resulting denominator singular
%Y = pzer(XY(rB+1:rB+cB,:),tolzero);   Deleted  09-Dec-2000
X = XY(1:rB,:); Y = XY(rB+1:rB+cB,:);
if issingular(Y,tol),
   error('axby0 system error: the resulting denominator singular.');
end;
X = X*(1/kappa);
X.v = Xv; Y.v = Xv;
X.h = Xh; Y.h = Xh;

%end .. axby0
