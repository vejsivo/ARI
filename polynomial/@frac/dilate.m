function B = dilate(A,k,h,var)
%DILATE    Dilate fraction
%
% Let A be a fraction in variable 'z^-1' or 'd'. The command
%     B = DILATE(A,K,H)
% where scalar integers K,H are dilating period, K>=1, and dilating
% phase, H>=0, returns fraction  B(z^-1) = z^-H * A(z^-K) .
% The arguments K,H are optional, its default value being K=2, H=0.
%
% For a fraction in variable 'z' or 'q', the numbering of H is
% reversed, the formula being  B(z) = z^-H * A(z^K) . This rule
% is for compatibility with TSP/DILATE.
%
% For a fraction in variable 's' or 'p', the dilatation is
% not defined.
%
% The variable symbol of B is taken from that of A. If A is constant
% and has no symbol then the symbol of B may be taken from an optional
% (fourth) argument VAR or from the global discrete time variable symbol.
%
% The meaning of the variable symbol has changed:
%  old var = new var to the K-th power
% The sampling period of B is K-times less than that of A. When VAR
% is applied, the sampling period of B is K-times less than that of VAR.
%
% See also POL/DILATE, TSP/DILATE.

%      Author: J. Jezek, 02-Jun-2000
%      Copyright(c) 2000 by Polyx, Ltd.
%      $ Revision $  $ Date 04-Oct-2000 $
%                    $ Date 14-Oct-2002 $
%                    $ Date 28-Feb-2003 $

global PGLOBAL;

ni = nargin;
if ni<4 |isempty(var), var = PGLOBAL.DISCRVAR;
end;
if ni<3 | isempty(h), h = 0;
end;
if ni<2 | isempty(k), k = 2;
end;
if ni<1,
   error('Not enough input arguments.');
end;

if ~isa(A,'frac'),
   error('Some argument but not 1st is invalidly fraction.');
end;
Avar = A.v;
if isempty(Avar), Avar = var;
end;
Bnum = 0; Bden = 0;
eval('Bnum = dilate(A.num,k,h,Avar); Bden = dilate(A.den,k,0,Avar);',...
   'error(peel(lasterr));');

if isa(A,'rdf'), B = Bnum/Bden;
elseif isa(A,'ldf'), B = Bden\Bnum;
elseif isa(A,'mdf'), B = Bnum./Bden;
elseif isa(A,'sdf'), B = Bnum*inv(Bden);
else
   error('Invalid 1st argument.');
end;

%end .. @frac/dilate
