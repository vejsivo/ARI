function gensym(varargin),
%GENSYM Set the global variable symbols for polynomial matrices
%
% GENSYM VALUE sets the value of the related global variable
% symbol to VALUE. Works with any number of input arguments such as
% GENSYM VALUE1 VALUE2 ... VALUEN. When called without arguments, 
% GENSYM sets the current values of the global variable symbol
% and continuous-time variable to 's' and  discrete-time variable to 'z'.
% It is the same as GENSYM S DISZ CONTS.
%
% The VALUE can be:
%    s  or  p                  Variable string for continuous-time,
%                                 derivative operator
%    z  or  q                  Variable string for discrete-time,
%                                 shift (advance) operator
%    z^-1, zi, or  d           Variable string for discrete-time,
%                                 shift (delay) operator
%    disz, discz, discrz       Default discrete-time variable
%    diszi, disczi, discrzi, 
%    disz^-1, discz^-1, 
%    discrz^-1
%    conts  or  contp          Default continuous-time variable
%
% See also GPROPS

%       Author(s):  S. Pejchova, M. Sebek 24-9-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 24-Sep-1998 15:35:34   $
%       $Revision: 3.0 $  $Date: 17-Jul-2000 16:45:34  S.Pejchova  $

global PGLOBAL;
na = nargin; test1 = 0;

eval('PGLOBAL.VARIABLE;', 'painit;');
if na,
  for ii=1:na,
    argm=varargin{ii};
    if isa(argm,'pol'),
       [as1,as2,ad]=size(argm);
       if all([as1,as2,ad]==1) & all(argm.c(:,:)==[0,1])
          PGLOBAL.VARIABLE = argm.v;
       else, test1=1;
       end;
    elseif isa(argm,'char') & ndims(argm)==2 & size(argm,1)==1,
      switch argm,
       case {'s','p','z^-1','d','z','q'}
          PGLOBAL.VARIABLE = argm;
       case 'zi',
          PGLOBAL.VARIABLE = 'z^-1';
       case {'disz','discz','discrz'}
          PGLOBAL.DISCRVAR = 'z';
       case {'diszi','disczi','discrzi','disz^-1','discz^-1','discrz^-1'}
          PGLOBAL.DISCRVAR = 'z^-1';
       case {'conts'}
          PGLOBAL.CONTVAR = 's';
       case {'contp'}
          PGLOBAL.CONTVAR = 'p';
       otherwise,
          test1=1;
       end;
    else,
       test1=1;
    end;
    if test1,
       error('Invalid input for global variable symbol.');
    end;
 end;
else,
   PGLOBAL.VARIABLE = 's';    % variable string
   PGLOBAL.DISCRVAR = 'z';    % discrete-time variable
   PGLOBAL.CONTVAR = 's';     % continuous-time variable
end;

%end .. gensym
