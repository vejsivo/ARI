function poldesk(argm1)
%POLDESK  Comprehensive hypertext documentation 
%and troubleshooting of the Polynomial Toolbox.
%
% POLDESK loads on-line versions of the Polynomial
% Toolbox 'Manual' and 'Commands' Help Desk Pages
% into Acrobat Reader.
% 
% Under Windows the macro looks for Acrobat Reader
% in the folder ..\polynomial\Pdf-files\Acrobat 4.0\Reader\
% that is provided as a part of the standard Polynomial 
% Toolbox installation. If the folder has been 
% deleted then the macro checks the file 
% C:\Program Files\Adobe\Acrobat 4.0\Reader\AcroRd32.exe 
% that typically results from an independent Acrobat 
% Reader installation. If also this file does not 
% exist then the user is asked to type the full path 
% name to the current Acrobat Reader location into
% a dialogue box.
%
% Under UNIX the macro calls the alias "acroread"
% that typically activates Acrobat Reader.
% Should such an alias not exist then it must be
% created by the user or by the system administrator.
%
% Typing POLDESK RECOVER or POLDESK('RECOVER') opens
% the dialogue box to type in the correct path name
% or alias.

%       Author(s):  S. Pejchova 30-3-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 23-Mar-2000 14:49:50   $
%                         $Date: 30-Aug-2002  S. Pejchova  $  

global PGLOBAL; msg1=0; msg2=0; aCRsTr=''; 

eval('PGLOBAL.FORMAT;', 'painit;');
str_pi=which('gprop.m');
fg1=find(str_pi==filesep);
if ~isempty(fg1), str_pi=str_pi(1:fg1(end));
else, str_pi='';
end;
str_xA=[str_pi,'private',filesep,'pacrpath'];
str_ver=version;
a1=str2num([str_ver(1),str_ver(3)]);
if a1 < 65,
    st1 = warning; warning off;
else, st1 = warning('query', 'all');
    warning off all;
end;
if nargin,
   if ~ischar(argm1),
      warning(st1);
      error('Invalid input argument.');
   elseif strcmp(argm1,'recover') 
      delete([str_xA,'.mat']); msg2=1;
   else, aCRsTr=argm1;
   end;
else,
   eval(['load(str_xA);'],'aCRsTr=''''; '); argm1='';  
end;

str_pdf=''; str_Rd='';
if (exist('Pdf-files','dir')==7),  str_pdf=[str_pi,'Pdf-files'];  end;

if (exist('OnLineManual.pdf','file')==2)&(exist('OnLineCommands.pdf','file')==2),
   str_OLM=which('OnLineManual.pdf');
   fg1=find(str_OLM==filesep);
   if ~isempty(fg1),
      str_pdf=str_OLM(1:fg1(end)); msg1=1;
   end; 
elseif ~isempty(str_pdf)
   addpath(str_pdf);
   if (exist('OnLineManual.pdf','file')==2)&(exist('OnLineCommands.pdf','file')==2),      
      str_OLM=[str_pdf,filesep,'OnLineManual.pdf']; msg1=1;
   end;
   if (exist('Acrobat 4.0','dir')==7),
      str_Rd=[str_pdf,filesep,'Acrobat 4.0',filesep,'Reader',filesep,'AcroRd32.exe'];
   end;
   rmpath(str_pdf);
end;

if ~msg1
   warning(st1);
   error(['Could not locate Polynomial Toolbox on-line PDF help files.' sprintf('\n') ...
          'Please  make  sure  that  the  files  OnLineManual.pdf and' sprintf('\n') ...
          'OnLineCommands.pdf  are  installed  correctly.']);
end;


if (isempty(aCRsTr))
   if isunix , 
      aCRsTr='acroread';    
   elseif isempty(str_Rd),
      aCRsTr=fullfile('C:','Program Files','Adobe',...
         'Acrobat 4.0','Reader','AcroRd32.exe');
   else, aCRsTr=str_Rd;
   end;
end
if ~isunix,
   if ~(exist(aCRsTr,'file')==2),
      fg1=find(aCRsTr==filesep); str_Us='';
      if ~isempty(fg1),  str_Us=aCRsTr(1:fg1(end)); end;
      if ~(exist(str_Us,'dir')==7),
         msg2=1;
      else,
         addpath(str_Us);
         if ~(exist(aCRsTr,'file')==2), msg2=1; end        
         rmpath(str_Us);
      end;
   end;
end;
warning(st1);
if msg2
   X_war2=figure('numbertitle','Off','Name','W A R N I N G',...
      'Units','normalized','MenuBar','none',...
      'Position',[0.24 0.6 0.7 0.28]);
%      'Position',[0.34 0.6 0.55 0.28],'WindowStyle','modal');
   uicontrol(X_war2,'Style','frame','Units','normalized',...
      'Position',[0.01 0.01 0.98 0.98]);
   
   war_str_n2=char('Could not locate an Acrobat Reader EXE file!',...
          ['Type the correct pathname, please.']);
   uicontrol(X_war2,'Style','text','Units','normalized',...
      'BackgroundColor',[0.5 0.55 0.7],...
      'ForegroundColor',[1 1 1],...
      'FontSize',10,'FontWeight','bold',...
      'Position',[0.1 0.55 0.8 0.32],'String',war_str_n2);
   uicontrol(X_war2,'Style','edit','Units','normalized',...
      'BackgroundColor',[1 1 1],...
      'ForegroundColor',[0 0 0],...
      'FontName','Courier','FontSize',9,...
      'Tag','Wr-Edx',...
      'Position',[0.1 0.3 0.8 0.24],'String',aCRsTr);
      
   uicontrol(X_war2,'Style','pushbutton','Units','normalized',...
      'BackgroundColor',[0.5 0.55 0.7],...
      'ForegroundColor',[1 1 1],...
      'String','OK','Position',[0.1 0.1 0.3 0.15],...
      'CallBack',...
      ['aCRsTrxX_X=get(findobj(''Tag'',''Wr-Edx''),''String''); ',...
      'close(gcf); pause(1); poldesk(aCRsTrxX_X); clear aCRsTrxX_X; ']);

   uicontrol(X_war2,'Style','pushbutton','Units','normalized',...
      'BackgroundColor',[0.5 0.55 0.7],...
      'ForegroundColor',[1 1 1],...
      'String','CANCEL','Position',[0.6 0.1 0.3 0.15],...
      'CallBack',['close(gcf);']);
else,
   fx1=find(aCRsTr==' ');  fx2=find(str_OLM==' ');
   if ~isempty(fx1), pli_xx=['"',aCRsTr,'"'];
   else, pli_xx=aCRsTr;
   end; 
   if ~isempty(fx2), str_OLM=['"',str_OLM,'"']; end;
   if nargin & ~strcmp(argm1,'recover'),
      eval(['save(str_xA,''aCRsTr'')';],'aCRsTr;'); 
   end;
   path(path);
   pause(1);
   eval(['! ',pli_xx,' ',str_OLM,' &'],'aCRsTr;'); pause(1);
end

%end .. poldesk
