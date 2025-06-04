function S = inv(R,arg2,arg3)
%INV    Inverse of matrix-den fraction
%
% For matrix-den fraction R, the command  S = INV(R)  returns
% the inverse of R. Matrix R must be square and nonsingular.
%
% An optional input argument TOL may specify the zeroing tolerance 
% to be used instead of the standard one.
%
% An optional input argument MET may specify the method used.
% The methods are  'def','int',  see POL/ADJ.

%         Author:  J. Jezek  07-Jan-2000
%         Copyright(c) by Polyx, Ltd.
%         $ Revision $  $ Date 26-Apr-2000 $
%                       $ Date 29-May-2000 $
%                       $ Date 06-Jul-2001 $
%                       $ Date 30-Sep-2002 $
%                       $ Date 14-Oct-2002 $

global PGLOBAL;

tol = PGLOBAL.ZEROING;
met = 'int';

if nargin>=2,
   if isa(arg2,'char'), met = arg2;
   elseif isa(arg2,'double'), tol = arg2;
   else error('Invalid 2nd argument.');
   end;
end;
if nargin==3,
   if isa(arg3,'char'), met = arg3;
   elseif isa(arg3,'double'), tol = arg3;
   else error('Invalid 3rd argument.');
   end;
end;

if ~isempty(tol),
   if length(tol)~=1 | ~isreal(tol) | tol<0 | tol>1,
      error('Invalid tolerance.');
   end;
end;   

Rs1 = R.frac.s(1); Rs2 = R.frac.s(2);
RRd = pol(zeros(Rs1,Rs1));
RRn = pol(zeros(Rs1,Rs2));
for i = 1:Rs1,
   [RRd(i,i),Delta] = plcm(R.frac.den(i,:),[],tol);
   RRn(i,:) = times(R.frac.num(i,:),Delta,tol);
end;
RRnad = 0; RRndet = 0;
eval('[RRnadj,RRndet] = adj(RRn,tol,met);', ...
   'error(peel(lasterr));');
if eq(RRndet,0,tol),
   error('Matrix is singular.');
end;
S = mdf(mtimes(RRnadj,RRd,tol),RRndet);

if strcmp(PGLOBAL.COPRIME,'cop'),
   S = coprime(S,tol);
end;
if strcmp(PGLOBAL.REDUCE,'red'),
   S = reduce(S);
else
   S = smreduce(S);
end;
if strcmp(PGLOBAL.DEFRACT,'defr'),
   S = defract(S);
end;

%end .. @mdf/inv
