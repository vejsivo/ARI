function sys = zpk(F,T,tol)
%ZPK   Convert two-sided polynomial to LTI
%     object in zero-pole-gain form
%
% For two-sided polynomial F, the command
%   YS = ZPK(F[,T[,TOL]])
% converts a system represented by F to a Control
% System Toolbox object in zero-pole-gain form.
%
% When T, the sampling time in the case of discrete-time
% systems, is not specified, it is taken from F.h . When
% F.h is empty, the default value 1 is taken. An optional
% tolerance may be specified in TOL. Its default value is
% the global zeroing tolerance.
%
% For backward conversion, see TSP/TSP.

%     Author:  J. Jezek, 22-Feb-2003
%     Copyright(c) 2000 by Polyx, Ltd.

if nargin<3, tol = [];
elseif ~isempty(tol),
   if ~isa(tol,'double'),
      error('Invalid tolerance.');
   end;
end;
if nargin<2, T = [];
elseif ~isempty(T),
   if ~isa(T,'double'),
      error('Invalid sampling period.');
      end;
end;

eval('sys = zpk(sdf(F),T,tol);', 'error(peel(lasterr));');

%end .. @tsp/zpk
