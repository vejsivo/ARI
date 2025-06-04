function pinit(varargin)
%PINIT  Initialize the Polynomial Toolbox
%
%   Modified version for ARI Course at CVUT FEL
%
% PINIT must be called at the beginning of every Polynomial 
% Toolbox session to set the global variable structure PGLOBAL
% to its default values.

%       Author(s):  S. Pejchova, M. Sebek 3-3-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 30-Apr-1999 11:46:50   $
%       $Revision: 3.0 $  $Date: 30-Jun-2000  J.Jezek  $
%                         $Date: 10-Jul-2000 15:08:00  S.Pejchova  $
%                         $Date: 02-Nov-2000  S.Pejchova  $
%                         $Date: 24-Jun-2001  J.Jezek  $
%                         $Date: 23-Sep-2001  J.Jezek  $
%                         $Date: 14-Mar-2002  Z.Hurak  $ 
%                         $Date: 27-Jan-2012  M.Sebek  $ modif. for ARI

pli_str=which('plic.mat');
if isempty(pli_str),error('Polynomial Toolbox file plic.mat does not exist.'); end;
try,
    load(pli_str),
catch,
    L_str='';
    L_num=[]; 
end
if (~exist('L_str','var')==1)|(~exist('L_num','var')), L_str=''; L_num=[]; end;
if isempty(L_num) | ~sum(L_num) | mod(sum(L_num),7),
   plicense('Init',varargin{1:nargin});
else,
   clear global PGLOBAL;
   global PGLOBAL;

   % Default values

   PGLOBAL.ZEROING = 1e-8;    % relative tolerance used for zeroing
   PGLOBAL.VERBOSE = 'no';    % flag to display extra comments during execution
   PGLOBAL.FORMAT = 'nice';   % display format
   PGLOBAL.ORDER = 'normal';  % display format order
   PGLOBAL.VARIABLE = 's';    % variable string
   PGLOBAL.DISCRVAR = 'z';    % discrete-time variable
   PGLOBAL.CONTVAR = 's';     % continuous-time variable
   PGLOBAL.REDUCE = 'nred';   % do not make the fraction reduced
                              %    after every operation
   PGLOBAL.COPRIME = 'ncop';  % do not make the fraction coprime
                              %    after every operation
   PGLOBAL.DEFRACT = 'ndefr'; % do not attempt to remove denominators
                              %    after every operation

   % New values
   na = nargin; test1 = 0;
   if na>=1,
      for ii=1:na,
         argm=varargin{ii};
         if isa(argm,'pol'),
            [vs1,vs2,vd] = size(argm);
            if all([vs1,vs2,vd]==1) & all(argm.c(:,:)==[0,1]),
            argm = argm.v;
            else
               error('Invalid input for global variable.');
            end;
         end;
         if isstr(argm) & ndims(argm)==2,
            switch argm,
            case {'s','p','z^-1','d','z','q'}
               PGLOBAL.VARIABLE = argm;
            case {'zi'}
               PGLOBAL.VARIABLE = 'z^-1';
            case {'symb', 'coef', 'rcoef', 'block','symbs','symbr','rootr','rootc','nice'}
               PGLOBAL.FORMAT = argm;
            case {'normal', 'reverse'}
               PGLOBAL.ORDER = argm;
            case {'no', 'yes'}
               PGLOBAL.VERBOSE = argm;
            case {'cop','coprime'}
               PGLOBAL.COPRIME = 'cop';
            case {'ncop','ncoprime'}
               PGLOBAL.COPRIME = 'ncop';
            case {'red','reduce'}
               PGLOBAL.REDUCE = 'red';
            case {'nred','nreduce'}
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
         elseif (isa(argm,'double')) & (all(size(argm)==1)) & ...
               isreal(argm) & argm>=0 & argm<=1,
            PGLOBAL.ZEROING = argm;
         else
            test1=1;
         end;
      end;
      if test1,
         error('Invalid input for global variable.');
      end;
   end; % for ii=1:na
   
   warning on;
   disp(' ');
   disp('  Polynomial Toolbox 3.0 initialized. To get started, type one of'); 
   disp('  these: helpwin  or poldesk.  For product information, visit');
   disp('  www.polyx.com  or  www.polyx.cz.');
   disp(' ');
   disp('  Special License for ARI Course at CVUT FEL');
   disp(' ');
end;

%end .. pinit


