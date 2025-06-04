function tolerance(tol),
%TOLERANCE  Set the global relative tolerance
%
% This commmand is used to set the global relative tolerance for 
% zeroing as follows:
%    TOLERANCE        Default. Same as TOLERANCE 1e-8
%    TOLERANCE 1e-8   Set global tolerance equal to 1e-8
%    TOLERANCE(1e-8)  Set global tolerance equal to 1e-8
%    TOLERANCE TOL    Set global tolerance equal to TOL
%    TOLERANCE(TOL)   Set global tolerance equal to TOL
% TOL is any real number


%       Author(s):  S. Pejchova, M. Sebek 29-5-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 29-May-1998 14:18:34   $
%       $Revision: 3.0 $  $Date: 18-Jul-2000 14:11:34  S.Pejchova  $
%                         $Date: 15-Jul-2001  J. Jezek  $
global PGLOBAL;
eval('PGLOBAL.ZEROING;', 'painit;');

if nargin == 0,
 PGLOBAL.ZEROING = 1e-8;
elseif ischar(tol) & ~isempty(tol),
   evalin('caller',['tolerance(',tol,');'],'error(lasterr)');
elseif (isa(tol,'double')) & (all(size(tol)==1)) & ...
      isreal(tol) & tol>=0 & tol<=1,
   PGLOBAL.ZEROING = tol;
else,
   error('Invalid input for global relative tolerance.');
end;

%end .. tolerance
