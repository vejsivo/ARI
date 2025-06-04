function painit(varargin)
%PAINIT  Initialize automatically the Polynomial Toolbox
%
% PAINIT is not supposed to be called by the user. It is called
% automatically when the user calls POL, TSP, RDF, LDF, MDF or SDF
% when the Polynomial Toolbox is not initialized. PAINIT sets
% the global variable structure PGLOBAL to its default values.
%
% See also: PINIT.

%       Author:  J. Jezek  14-Apr-2000
%       Copyright (c) 2000 by Polyx, Ltd.
%       $ Revision $  $ Date 24-Jun-2001   J.Jezek $
%                     $ Date 23-Sep-2001   J.Jezek $

pli_str=which('plic.mat');
if isempty(pli_str),error('Polynomial Toolbox file plic.mat does not exist.'); end;

eval(['load ',pli_str],'L_str='''';, L_num=[];');
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
   PGLOBAL.ORDER = 'normal';  % display format orde   
   PGLOBAL.VARIABLE = 's';    % variable string
   PGLOBAL.DISCRVAR = 'z';    % discrete-time variable
   PGLOBAL.CONTVAR = 's';     % continuous-time variable
   PGLOBAL.REDUCE = 'nred';   % do not make the fraction reduced
                              %    after every operation
   PGLOBAL.COPRIME = 'ncop';  % do not make the fraction coprime
                              %    after every operation
   PGLOBAL.DEFRACT = 'ndefr'; % do not attempt to remove denominators                              
                              %    after every operation
                              
   warning on;                           
   disp(' ');
   disp('  Polynomial Toolbox initialized.');
   disp(' ');
end;

%end .. painit


