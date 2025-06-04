function t = checkpb
%CHECKPB  Check for conflicts of Polynomial Toolbox variables and existing variables.
%
% CHECKPB returns:
%   T = [] if the Polynomial Toolbox is set correctly
%   T = 1  if the Polynomial Toolbox is not initialized or the
%          global property structure is incorrect
%   T = 2  if the global property VARIABLE has an incorrect value
%   T = 3  if the global property TOLERANCE has an incorrect value
%   T = 4  if the global property VERBOSE has an incorrect value
%   T = 5  if the global property FORMAT has an incorrect value
%   T = 6  if there already exists a variable s, p, zi, d, z, q or v
%          and the polynomial function of the same name is not active
%   T = 7  if the global property DISCRVAR has an incorrect value
%   T = 8  if the global property CONTVAR has an incorrect value
%   T = 9  if the global property REDUCE has an incorrect value
%   T = 10 if the global property COPRIME has an incorrect value
%   T = 11 if the global property DEFRACT has an incorrect value
% If there are more problems at the same time then the macro returns
% a vector listing the 'error codes' above.
%
% If called without output argument then CHECKPB displays a status message.

  
%       Author(s):  S. Pejchova 09-11-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 20-Nov-1998 12:42:50   $
%       $Revision: 3.0 $  $Date: 20-Jul-2000 15:14:34  S.Pejchova  $
 

global PGLOBAL; war_str=''; test1=0; test2=[];
eval('PGLOBAL.VARIABLE;','test1=1;');
if ~test1,
   names=fieldnames(PGLOBAL);
   for ii=1:length(names),
      switch names{ii},
       case 'VARIABLE',
         argm=PGLOBAL.VARIABLE;
         if ~ischar(argm)|...
            isempty(strmatch(argm,{'s';'p';'z^-1';'d';'z';'q'},'exact')),
            war_str=char(war_str,'*Global property VARIABLE has incorrect value.',...
               ' Use PINIT, GENSYM or GPROPS to recover.');
            test2=[test2,2];
         end;
       case 'ZEROING',
         argm=PGLOBAL.ZEROING;
         if ~isa(argm,'double') | (any(size(argm)>1)),
            war_str=char(war_str,'*Global property TOLERANCE has incorrect value.',...
               ' Use PINIT, TOLERANCE or GPROPS to recover.');
             test2=[test2,3];
        end;
       case 'VERBOSE',
         argm=PGLOBAL.VERBOSE;
         if ~ischar(argm)|...
            isempty(strmatch(argm,{'no';'yes'},'exact')),
            war_str=char(war_str,'*Global property VERBOSE has incorrect value.',...
                   ' Use PINIT, VERBOSE or GPROPS to recover.');
            test2=[test2,4];
         end;
       case 'FORMAT',
         argm=PGLOBAL.FORMAT;
         if ~ischar(argm)|...
            isempty(strmatch(argm,{'symb';'coef';'rcoef';'block';'symbs';'symbr';'rootr';'rootc'},'exact')),
            war_str=char(war_str,'*Global property FORMAT has incorrect value.',...
                   ' Use PINIT, PFORMAT or GPROPS to recover.');
            test2=[test2,5];
         end;
       case 'DISCRVAR',
         argm=PGLOBAL.DISCRVAR;
         if ~ischar(argm)|...
            isempty(strmatch(argm,{'z';'z^-1'},'exact')),
            war_str=char(war_str,'*Global property DISCRVAR has incorrect value.',...
                   ' Use PINIT, GENSYM or GPROPS to recover.');
            test2=[test2,7];
         end;
       case 'CONTVAR',
         argm=PGLOBAL.CONTVAR;
         if ~ischar(argm)|...
            isempty(strmatch(argm,{'s';'p'},'exact')),
            war_str=char(war_str,'*Global property CONTVAR has incorrect value.',...
                   ' Use PINIT, GENSYM or GPROPS to recover.');
            test2=[test2,8];
         end;
       case 'REDUCE',
         argm=PGLOBAL.REDUCE;
         if ~ischar(argm)|...
            isempty(strmatch(argm,{'red';'nred'},'exact')),
            war_str=char(war_str,'*Global property REDUCE has incorrect value.',...
                   ' Use PINIT, or GPROPS to recover.');
            test2=[test2,9];
         end;
       case 'COPRIME',
         argm=PGLOBAL.COPRIME;
         if ~ischar(argm)|...
            isempty(strmatch(argm,{'cop';'ncop'},'exact')),
            war_str=char(war_str,'*Global property COPRIME has incorrect value.',...
                   ' Use PINIT, or GPROPS to recover.');
            test2=[test2,10];
         end;
       case 'DEFRACT',
         argm=PGLOBAL.DEFRACT;
         if ~ischar(argm)|...
            isempty(strmatch(argm,{'defr';'ndefr'},'exact')),
            war_str=char(war_str,'*Global property DEFRACT has incorrect value.',...
                   ' Use PINIT, or GPROPS to recover.');
            test2=[test2,11];
         end;
        otherwise, test1=1;
      end;
   end;
end;
if test1,
   war_str=char(war_str,'*Global property structure incorrect. Use PINIT!');
   test2=[1,test2];
end;
save vAr_TmxX_b war_str test1;
evalin('caller','save vAr_TmxX_a;');
load vAr_TmxX_a; Sym_ar={'s';'p';'zi';'d';'z';'q';'v'};
war_str0='';  war_str1=['*Polynomial function '];
war_str2=[' is not active.'];
war_str3=[' Clear or rename the variable '];
war_str4=[' to activate the function.'];
for ii=1:7,
   if (exist(Sym_ar{ii})==1),
      war_str0=char(war_str0,[war_str1,Sym_ar{ii},war_str2],...
         [war_str3,Sym_ar{ii},war_str4]);
   end;
end;
load vAr_TmxX_b;
delete vAr_TmxX_a.mat; delete vAr_TmxX_b.mat;
if ~isempty(war_str0),
   war_str0=war_str0([2:end],:);
   war_str=char(war_str,war_str0,''); 
   test2=[test2,6];
elseif ~isempty(war_str),
   war_str=char(war_str,'');
else,
   war_str=char('',['*Polynomial Toolbox is correctly set.'],'');
end;
if nargout,
   t=test2;
else,
   disp(war_str);
end;

%end .. checkpb
