function [L_str,L_num]=plicense(argm1,varargin)
%PLICENSE  Polynomial Toolbox license number.
%
%   PLICENSE returns a string containing the Polynomial
%   Toolbox license number
  
%       Author(s):  S. Pejchova 19-04-99
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 03-June-1999 11:03:50   $
%       $Revision: 3.0 $  $Date: 18-Jul-2000 14:17:34  S.Pejchova  $

pli_str=which('plic.mat');
if isempty(pli_str),error('Polynomial Toolbox file plic.mat does not exist.'); end;

if nargin,
   if ~ischar(argm1), error('Invalid argument.'); end;
   if ~strcmp(argm1,'Init'), error('Invalid argument.'); end;
   
   lic_str_n1=['\n  Type your Polynomial Toolbox personal license number, please: \n \n  '];
   lic_str_n2=['\n  Invalid Polynomial Toolbox personal license number.\n',...
               '  Try once more or type Q for quit. \n \n  '];
   Lic_str=lic_str_n1; msg_err=0;
   while ~msg_err,
      L_str=input(Lic_str,'s'); L_num=[];
      if isempty(L_str),
        pause(1);
        Lic_str=lic_str_n2; msg_err=0;          
      elseif strcmpi(L_str(1),'q'), msg_err=2;
      else,
        L_str=strrep(L_str,' ',''); str_l1=strrep(L_str,'-','');
        if (length(str_l1)>9)&(length(str_l1)<16),
           eval(['L_num=str2num((str_l1)'')''; msg_err=1;'],'msg_err=0;'),
        end;
        L_sum=sum(L_num);   
        if isempty(L_sum)|(L_sum==0)|(mod(L_sum,7)), msg_err=0; end;
        pause(1);
        Lic_str=lic_str_n2; 
      end;
   end;
      
   if msg_err==1,
      str_at2=['error(''Permission to write in the Polynomial Toolbox denied.'');'];
      fx=find(pli_str==' ');
      if ~isempty(fx), pli_xx=['"',pli_str,'"'];
      else, pli_xx=pli_str;
      end;
      if isunix,
         eval(['! chmod +w ',pli_xx],str_at2);
      else,
         eval(['! attrib -r ',pli_xx],str_at2);
      end;     
      str_at1=['save(pli_str,''L_str'',''L_num'');'];
      eval(str_at1,str_at2);
      
      if nargin==1,  pinit;
      else,
         pinit(varargin{1:nargin-1});
      end;
      
      % simulink
      eval('vs = ver(''simulink'');', 'vs.Version = ''2.0'';');
      vsVer = str2num(vs.Version);
      
      if vsVer >= 3,
         disp(sprintf('\n'));
         warning(sprintf([' You seem to be using the version 3 of Simulink.\n', ...
               '  If it is true (you can check it by typing VER SIMULINK),\n', ... 
               '  replace the file\n\n', ...
            	'	        .../polynomial/polblock.mdl\n\n', ... 
            	'  in the Polynomial Toolbox 2.0 main directory \n', ...
            	'  by the file\n\n  ', ... 
         		'	        .../polynomial/simulink3/polblock.mdl \n\n', ...
      			'  if you have not done it yet.\n']));
      end;   
      % ^^ simulink ^^

   end; 
else,
   L_str=''; L_num=[];
   eval(['load(pli_str);'],'L_str='''';, L_num=[];');
end;  
   
%end .. plicense
