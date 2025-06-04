% ----------------------------------------------------------------------------
% POLYNOMIAL TOOLBOX 2.0.0 RELEASE NOTES                  27-Apr-99
% ----------------------------------------------------------------------------
% CONTENTS
% * Installation instructions
% * Documentation
% * A note for SIMULINK 3 users on Win95/98/NT
% * PolyX Ltd.
% ----------------------------------------------------------------------------
% INSTALLATION INSTRUCTIONS
% 
% The Polynomial Toolbox can be installed in the following simple steps:
% 
% Win 95/98/NT users:
% * Copy the whole folder  \polynomial  including all the subfolders from the 
%   Polynomial Toolbox CD-ROM to your PC, preferably next to other MATLAB 
%   toolboxes that all are placed in the   …\MATLAB\toolbox  or
%    …\MATLABR11\toolbox folders.
% * Add the folder \polynomial to your MATLAB path (for instance by using the
%   MATLAB Path Browser).
% * If you use version 3 of SIMULINK, replace the file 
%   ...\polynomial\polblock.mdl, residing in the main Polynomial Toolbox 
%   directory, by the file  ...\polynomial\simulink3\polblock.mdl. 
%   The current version of SIMULINK can be checked by typing "ver simulink" 
%   in the MATLAB main window. 
% * You are recommended to add to your  startup.m  file a new line 
%   containing the command PINIT. With this modification the Polynomial 
%   Toolbox is automatically initialized at the beginning of every MATLAB 
%   session. If you do not do this then you will have to type PINIT manually 
%   each time you start a Polynomial Toolbox session.
% * When using the Polynomial Toolbox for the first time after the installation
%    you will be asked to provide your personal license number.
% * The standard configuration of the Polynomial Toolbox contains an Acrobat
%   Reader placed in  the folder \polynomial\Pdf-files\Reader, which guarantees
%   easy use of the on-line documentation by the POLDESK command. This 
%   configuration requires no action during the installation and it is 
%   recommended for most users. For more details on using the on-line 
%   documentation see the section DOCUMENTATION below.
% 
% UNIX users:
% * Copy the whole directory /polynomial  including all subdirectories from 
%   the Polynomial Toolbox CD-ROM to your system, preferably next to the other 
%   MATLAB toolboxes that all are placed in  …/MATLAB/toolbox   or  
%    …/MATLABR11/toolbox  directories.
% * Add the directory /polynomial to your MATLAB path.
% * If you use version 3 of SIMULINK, replace the file 
%   .../polynomial/polblock.mdl, residing in the main Polynomial Toolbox 
%   directory, by the file  .../polynomial/simulink3/polblock.mdl. 
%   The current version of SIMULINK can be checked by typing "ver simulink" 
%   in the MATLAB main window. 
% * You are recommended to add to your startup.m file a new line containing 
%   the command PINIT. With this modification the Polynomial Toolbox is 
%   automatically initialized at the beginning of every MATLAB session. If you 
%   do not do this then you will have to type PINIT manually each  time you 
%   start a new Polynomial Toolbox session.
% * When using the Polynomial Toolbox for the first time after the installation
%   you will be asked to provide your personal license number.
% * To access the Polynomial Toolbox on-line documentation by the command
%   POLDESK your UNIX system is supposed to run Acrobat Reader by the
%   usual command ‘acroread’. If this is not the case then you must create such 
%   an alias, or ask your system administrator for help. You may wish to save 
%   disk space by deleting the directory /polynomial/Pdf-files/Reader that is 
%   of no use for UNIX systems that have a standard installation of Acrobat
%   Reader. For more details on using the on-line documentation see the section
%   DOCUMENTATION below. 
% 
% -----------------------------------------------------------------------------
% DOCUMENTATION
% 
% There are two document volumes provided with the Polynomial Toolbox:  
% MANUAL and COMMANDS. The printed MANUAL volume is delivered 
% with the Polynomial Toolbox CD-ROM. The printed COMMANDS volume 
% may be purchased separately (see the PolyX Web site or contact 
% info@polyx.cz). 
% 
% Ready-to-print electronic versions of the both MANUAL and COMMANDS 
% are also available. They may be found in  \polynomial\Pdf-files in the files 
% manual.pdf and commands.pdf , which are readable by Acrobat Reader. Users 
% are allowed to print these files for their own use but should not distribute 
% them any further. For more copyright details see the License Agreement. 
% 
% On-line electronic versions of the MANUAL and COMMANDS are also provided. 
% They are located in \polynomial\Pdf-files in the files  OnLineManual.pdf 
% and  OnLineCommands.pdf . They are normally accessed by the Polynomial 
% Toolbox command POLDESK but users are free to create other arrangements.
% 
% In Win95\98\NT, POLDESK by default uses Acrobat Reader located in the 
% Polynomial Toolbox  folder \polynomial\Pdf-files\Reader . This configuration
% is generally recommended. If, however, an experienced user wishes to employ 
% a different version of Acrobat Reader located elsewhere then the entire folder 
% \polynomial\Pdf-files\Reader may simply be deleted. During the next execution
% POLDESK will look for Acrobat Reader in the standard location 
% C:\Acrobat3\Reader\AcroRd32.exe or will ask the user to provide a valid path
% name.
% 
% In UNIX, first of all, the whole directory  /polynomial/Pdf-files/Reader is 
% useless and can be deleted if a standard installation of Acrobat Reader 
% exists. POLDESK by default calls the command ‘acroread’ that typically runs 
% Acrobat Reader on a UNIX system. If this is not the case the user or a system 
% administrator can create such an alias.
% 
% Alternatively, the user of each system can type POLDESK RECOVER that opens
% a dialogue window where the user can type in a valid pathname.
% ------------------------------------------------------------------------------
% A NOTE FOR SIMULINK 3 USERS ON WIN95/98/NT
% 
% Under the MS Windows operating systems the way the "simulink" 
% command is processed slightly differs in versions 2 and 3 of 
% SIMULINK. The instructions in the Polynomial Toolbox 2.0 Manual 
% (pages 83-84) refer to version 2 of SIMULINK. If you use SIMULINK 3 
% under Windows then please proceed in one of the two following ways:
% 
% EITHER:
% 
% type 
% 
% >>simulink 
% 
% to open the Simulink Library browser. The Polynomial Toolbox 2.0  
% Simulink library now is directly accessible within the browser 
% along with the other Simulink libraries.
% 
% OR:
% 
% type 
% 
% >>simulink3
% 
% to open the Simulink Library window. Follow the instructions 
% in the Polynomial Toolbox 2.0 Manual, pages 83-84.
% 
% For further information consult the SIMULINK manual 
% (Using Simulink, Version 3).
% ------------------------------------------------------------------------------
% POLYX, LTD
% ----------
% ADDRESS: Jarni 4, 16000 Praha 6, Czech Republic
% PHONE  +420-603-844561,    FAX    +4202-6884554 
% WEBS:   www.polyx.cz       OR     www.polyx.com
% 
% E-MAIL BOXES:
% info@polyx.cz   	general questions and suggestions
% sales@polyx.cz  	purchase related questions 
% support@polyx.cz 	technical questions
% bugs@polyx.cz 	        bug reports
% doc@polyx.cz 	        documentation related questions and error reports
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
% Copyright (c) 1999 PolyX, Ltd. All Rights Reserved.
% Version 2.0.0  27-Apr-99
% ------------------------------------------------------------------------------
