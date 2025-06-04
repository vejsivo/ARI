function gprops(varargin)
%GPROPS   Display/modify global properties
%
% GPROPS VALUE sets the value of the related global property
% to VALUE. Works with any number of input arguments such as
% GPROPS VALUE1 VALUE2 ... VALUEN. When called without arguments, 
% GPROPS displays the current values of all the global properties 
% along with the list of alternatives. 
  
%       Author(s):  S. Pejchova, M. Sebek 3-3-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 10-Nov-1998 14:17:50   $
%       $Revision: 3.0 $  $Date: 06-Apr-2000  J.Jezek   $ 
%                         $Date: 10-Jul-2000  S.Pejchova   $ 
%                         $Date: 02-Nov-2000  S.Pejchova   $ 
%                         $Date: 23-Sep-2001  J.Jezek   

global PGLOBAL;
eval('PGLOBAL.VARIABLE;', 'painit;');

na = nargin; test1 = 0;
if na
  for ii=1:na,
    argm=varargin{ii};
    if isa(argm,'pol'),
       [as1,as2,ad]=size(argm);
       if all([as1,as2,ad]==1) & all(argm.c(:,:)==[0,1])
          PGLOBAL.VARIABLE = argm.v;
       else test1=1;
       end;
    elseif isa(argm,'char') & ndims(argm)==2 & size(argm,1)==1,
      switch argm,
       case {'s','p','z^-1','d','z','q'}
          PGLOBAL.VARIABLE = argm;
       case 'zi',
          PGLOBAL.VARIABLE = 'z^-1';
       case {'symb', 'coef', 'rcoef','block','symbs','symbr','rootr','rootc','nice'}
          PGLOBAL.FORMAT = argm;
       case {'normal', 'reverse'}
          PGLOBAL.ORDER = argm;
       case {'no', 'yes'}
          PGLOBAL.VERBOSE = argm;
       case {'cop','coprime'}
          PGLOBAL.COPRIME = 'cop';
       case {'ncop','ncoprime'}
          PGLOBAL.COPRIME = 'ncop';
       case {'red','reduce','reduced'}
          PGLOBAL.REDUCE = 'red';
       case {'nred','nreduce','nreduced'}
          PGLOBAL.REDUCE = 'nred';
       case {'def','defr','defrac','defract'}
          PGLOBAL.DEFRACT = 'defr';
       case {'ndef','ndefr','ndefrac','ndefract'}
          PGLOBAL.DEFRACT = 'ndefr';
       case {'disz','discz','discrz'}
          PGLOBAL.DISCRVAR = 'z';
       case {'diszi','disczi','discrzi','disz^-1','discz^-1','discrz^-1'}
          PGLOBAL.DISCRVAR = 'z^-1';
       case {'conts','contp'}
          PGLOBAL.CONTVAR = argm(5);
       otherwise
          argm = str2num(argm);
          if ~isempty(argm) & length(argm)==1 & ...
                isreal(argm) & argm>=0 & argm<=1,
             PGLOBAL.ZEROING = argm;
          else test1 = 1;
          end;
      end;  %switch argm
    elseif (isa(argm,'double')) & length(argm)==1  & ...
            isreal(argm) & argm>=0 & argm<=1,
       PGLOBAL.ZEROING = argm;
    else test1=1;
    end;
    if test1,
       error('Invalid input for global variable.');
    end;
 end; % for ii=1:na
 
else, % no arguments
  disp(' ');
  disp('Global polynomial properties:');
  disp(' ');
  cop = PGLOBAL.COPRIME;
  red = PGLOBAL.REDUCE;
  defr = PGLOBAL.DEFRACT;
  if strcmp(PGLOBAL.DISCRVAR,'z'), dvar = 'disz';
  else dvar = 'disz^-1';
  end;
  if strcmp(PGLOBAL.CONTVAR,'s'), cvar = 'conts';
  else cvar = 'contp';
  end;
  
  col1=char('PROPERTY NAME:','variable symbol','zeroing tolerance   ',...
     'verbose level','display format',' ','display order',...
     'coprime flag','reduce flag','defract flag',...
     'discrete variable','continuous variable');
  col2=char('CURRENT VALUE:  ',PGLOBAL.VARIABLE,num2str(PGLOBAL.ZEROING),...
     PGLOBAL.VERBOSE,PGLOBAL.FORMAT,' ',PGLOBAL.ORDER,...
     cop,red,defr,dvar,cvar);
  col3=char('  AVAILABLE VALUES:',[' ''s'',''p'',''z^-1'',''d'',''z'',''q'''],...
            '  any real number from 0 to 1',[' ''no'', ''yes'''],...
            [' ''symb'', ''symbs'', ''nice'', ''symbr'', ''block'''],...
            [' ''coef'', ''rcoef'', ''rootr'', ''rootc'''],...
            [' ''normal'', ''reverse'''],...
            [' ''cop'', ''ncop'''],...
            [' ''red'', ''nred'''],...
            [' ''defr'', ''ndefr'''],...
            [' ''disz'', ''disz^-1'''],...
            [' ''conts'', ''contp''']);
  disp([col1,col2,col3]);
  disp(' ');
end;

%end .. gprops
