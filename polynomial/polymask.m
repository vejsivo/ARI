function [A,B,C,D,Label,var] = polymask(P,Q,Type,FVar, NStat);
% POLYMASK  Initialize the mask of Polynomial Library for Simulink.
% 
% This function is intended to be called only by the POL Block in Simulink.

%	Author(s): M. Hromcik, M. Sebek 7-9-98
%	Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 7-Sep-1998 10:28:34   $

% Effect on other properties:
% A,B,C,D are standard Matlab matrices.

% variable assigning
Pv = P.var;
Qv = Q.var;

switch FVar,		% forcing variable
  case 1,		% no force - unify variables
   if ~strcmp(Qv, Pv),
      if isempty(Qv),
       Q.var = Pv;
       var = Pv;
      elseif isempty(Pv),
       P.var = Qv;   
       var = Qv;
      elseif (strcmp(Qv, 's') & strcmp(Pv, 'p')) | ...	   % 's', 'p' => 's'
             (strcmp(Pv, 's') & strcmp(Qv, 'p')),
       P.var = 's';
       Q.var = 's';
       var = 's';
      elseif (strcmp(Qv, 'z') & strcmp(Pv, 'q')) | ...     % 'z', 'q' => 'z'
             (strcmp(Pv, 'z') & strcmp(Qv, 'q')),     
       P.var = 'z';
       Q.var = 'z';  
       var = 'z';
      elseif (strcmp(Qv, 'z^-1') & strcmp(Pv, 'd')) | ...  % 'd', 'z^-1' => 'd'
             (strcmp(Pv, 'z^-1') & strcmp(Qv, 'd')),     
       P.var = 'd';
       Q.var = 'd'; 
       var = 'd';
      
      else 		   % conflict
       var = '';
     end;
     
   else
     var = Pv; 
   end;
  case 2,		% force 's'
   var = 's';
   P.var = 's';
   Q.var = 's';
  case 3,		% force 'p'
   var = 'p';
   P.var = 'p';
   Q.var = 'p';
  case 4,		% force 'z'
   var = 'z';
   P.var = 'z';
   Q.var = 'z';
  case 5,		% force 'q'
   var = 'q';
   P.var = 'q';
   Q.var = 'q';
  case 6,		% force 'd'
   var = 'd';
   P.var = 'd';
   Q.var = 'd';
  case 7,		% force 'z^-1'
   var = 'z^-1';
   P.var = 'z^-1';
   Q.var = 'z^-1';
  otherwise 
   error('Unknown command option.');
end;	%switch 

if ~isempty(var),
  var = ['(' var ')'];
end;  

% d-to-z conversion - built in "...2s" functions.
%if strcmp(var,'(d)') | strcmp(var,'(z^-1)'), 
%  warning off;
%  [P,Q] = dz(P,Q);
%  warning on;
%end;

% mf-to-state-space conversion + mask string creating 

if Type == 1,
  %eval('[A, B, C, D] = lmf2ss(P,Q);', 'errordlg(lasterr); error(lasterr);');
  [A, B, C, D] = lmf2ss(P,Q);
  opt = ' Q\P ';
  SumDeg = sum(deg(Q,'row'));

elseif Type == 2,
  [A, B, C, D] = rmf2ss(P,Q);
  opt = ' P/Q ';
  SumDeg = sum(deg(Q,'col'));
 
else 
  error('Unknown command option.');
end;

if strcmp(var, '(z^-1)'),
  opt = opt(2:4);
end;
Label = ['y' var ' =' opt 'u' var];

% Check # initial states:
if ~(SumDeg == NStat  |  NStat == 0),
  error(['Invalid setting - the number of initial states must equal ',...
        num2str(SumDeg)]);
end;

%end ../polymask  

