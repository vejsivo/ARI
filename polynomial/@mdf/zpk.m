function sys = zpk(F,T,tol)
%ZPK   Convert matrix-den fraction to LTI object
%      in zero-pole-gain form
%
% For matrix-den fraction F, the command
%    SYS = ZPK(F[,T,[TOL]])
% converts a system represented by the fraction to a Control
% System Toolbox object in zero-pole-gain form.
%
% When T, the sampling time in the case of discrete-time
% systems, is not specified, it is taken from F.h . When
% F.h is empty, the default value 1 is taken. An optional
% tolerance may be specified in TOL. Its default value is
% the global zeroing tolerance.
%
% For backward conversion, see MDF/MDF.

%     Author:  J. Jezek, 16-Feb-2000
%     Copyright(c) 2000 by Polyx, Ltd.
%     $ Revision $  $ Date 26-Apr-2000 $
%                   $ Date 31-May-2000 $
%                   $ Date 14-Oct-2002 $

global PGLOBAL;

if nargin>=2,
   if ~isa(T,'double') | length(T)~=1 | ~isreal(T) | T<0,
      error('Invalid sampling period.');
   end;
   if nargin==3,
      if ~isa(tol,'double') | length(tol)~=1 | ...
            ~isreal(tol) | tol<0 | tol>1,
         error('Invalid tolerance.');
      end;
   else
      tol = PGLOBAL.ZEROING;
   end;
else
   if ~isempty(F.frac.h), T = F.frac.h;
   else T = 1;
   end;
   tol = PGLOBAL.ZEROING;
end;    % TOL is accepted but not used 

N = F.frac.num; D = F.frac.den;
Nd = N.d; Dd = D.d;
Fs1 = F.frac.s(1); Fs2 = F.frac.s(2);
num = cell(Fs1,Fs2); den = cell(Fs1,Fs2);
for i = 1: Fs1,
   for j = 1:Fs2,
      num{i,j} = squeeze(N.c(i,j,Nd+1:-1:1)).';
      den{i,j} = squeeze(D.c(i,j,Dd+1:-1:1)).';
   end;
end;

[Z,KZ] = pol2root(mat2pol(num));
[P,KP] = pol2root(mat2pol(den));
K = KZ ./ KP;

Imp = 'Control Toolbox not found.';
switch F.frac.v
case {'s','p'}
   eval('sys = zpk(Z,P,K);','disp(Imp)');
otherwise
   eval('sys = zpk(Z,P,K,T);','disp(Imp)');
end;

%end .. @mdf/zpk
