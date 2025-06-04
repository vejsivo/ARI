function pver
%PVER  Polynomial Toolbox version information
%
% The commmand   PVER 
% displays the current Polynomial Toolbox version number.
  
%       Author(s):  S. Pejchova 19-11-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 29-Apr-1999 08:36:50   $
%       $Revision: 3.0 $  $Date: 18-Jul-2000 12:00:34  S.Pejchova  $
 

global PGLOBAL;
eval('PGLOBAL.VARIABLE;','painit;');
disp('---------------------------------------------------');
disp(['POLYNOMIAL TOOLBOX  Version ',pversion,' on ',computer]);
if (exist('plicense','file')==2) | (exist('plicense','file')==6),
   disp(['POLYNOMIAL TOOLBOX  License Number: ',plicense]);
end;
disp('---------------------------------------------------');

%end .. pver
