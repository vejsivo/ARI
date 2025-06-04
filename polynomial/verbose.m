function verbose(string),
%VERBOSE  Set the global verbose level in the Polynomial Toolbox
%
% This commmand is used to switch between different global verbose 
% levels as follows:
%    VERBOSE        Default. Same as VERBOSE NO
%    VERBOSE NO     No comments during the execution
%    VRBOSE YES     First level of comments during execution

%       Author(s):  S. Pejchova, M. Sebek 29-5-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 25-Sep-1998 10:00:34   $
%       $Revision: 3.0 $  $Date: 18-Jul-2000 12:05:34  S.Pejchova  $

global PGLOBAL;
eval('PGLOBAL.VERBOSE;', 'painit;');

if nargin,
   if ischar(string) & (strcmp(lower(string),'no')| strcmp(lower(string),'yes')),
      PGLOBAL.VERBOSE = lower(string);
   else,
      error('Invalid verbose level.');
   end; 
else,
   PGLOBAL.VERBOSE = 'no';
end;

%end .. verbose
