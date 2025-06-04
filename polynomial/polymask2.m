function [A,B,C,D,Label,var] = polymask(Fstr,FVar,NStat);
% POLYMASK  Initialize the mask of Polynomial Library for Simulink.
% 
% This function is intended to be called only by the POL Block in Simulink.

%	Author(s): M. Hromcik 6-12-2000
%	Copyright (c) 1998 by Polyx, Ltd.
%  Version 3

F = evalin('base',Fstr);

% variable assigning
Fv = F.var;

switch FVar,		% forcing variable
case 1,		% no force - leave variable unchanged
   var = Fv;
case 2,		% force 's'
   var = 's';
   F.var = 's';
case 3,		% force 'p'
   var = 'p';
   F.var = 'p';
case 4,		% force 'z'
   var = 'z';
   F.var = 'z';
case 5,		% force 'q'
   var = 'q';
   F.var = 'q';
case 6,		% force 'd'
   var = 'd';
   F.var = 'd';
case 7,		% force 'z^-1'
   var = 'z^-1';
   F.var = 'z^-1';
otherwise 
   error('Invalid command option.');
end;	%switch 

if ~isempty(var),
  var = ['(' var ')'];
end;  

% mf-to-state-space conversion + mask string creating 

[A,B,C,D] = abcd(F);
SumDeg = size(A,1);

if length(Fstr)<10,
   Label = Fstr;
   if length(Fstr)<=4,
      Label = ['y' var ' = ' Label ' u' var];
   end;
   
else
	if isa(F,'ldf'),
    	opt = ' D\N ';
   elseif isa(F,'rdf'),
   	opt = ' N/D ';
	elseif isa(F,'mdf'),
   	opt = ' N./D ';
	elseif isa(F,'sdf'),
   	opt = ' N/d ';
	end;
   Label = ['y' var ' =' opt 'u' var];
end;

% Check # initial states:
if ~(SumDeg == NStat  |  NStat == 0),
  error(['Invalid setting - the number of initial states must equal ',...
        num2str(SumDeg)]);
end;

%end .. polymask
