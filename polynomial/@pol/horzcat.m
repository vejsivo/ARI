function C = horzcat(varargin)
%HORZCAT   Horizontal concatenate polynomials
%
% The commmand
%    C = [A B] 
% returns the horizontal concatenation of the polynomial matrices 
% A and B.  A and B must have the same numbers of rows.  
%
% Any number of matrices may be concatenated within 
% one pair of brackets. Horizontal and vertical concatenation may be 
% combined together as for standard MATLAB matrices. The commmand
%    C = [A,B]  or  C = HORZCAT(A,B)
% works similarly.
%  
%  See also POL/VERTCAT.

%       Author(s): M. Hromcik, M. Sebek, S. Pejchova 16-2-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 21-Apr-1998 10:44:34   $
%       $Revision: 3.0 $  $Date: 11-Aug-1999 12:00:00   J. Jezek  $
%                         $Date: 11-Oct-1999 12:00:00   J. Jezek  $
%                         $Date: 16-Jun-2000 12:00:00   J. Jezek  $
%                         $Date: 03-Aug-2000 09:30:00   J. Jezek  $
%                         $Date: 07-Nov-2000 12:00:00   J. Jezek  $
%                         $Date: 25-Jan-2002            J. Jezek  $
%                         $Date: 28-Feb-2003            J. Jezek  $

% Effect on other properties:
% C.v: Variable is kept if it is the same in all input matrices, 
%      otherwise it is set equal to the current value of PGLOBAL.VARIABLE.
%      However, if there are both 'z' and 'z^-1' among the input variables
%      (and no others), the result is a two-sided polynomial object.
% C.u: UserData are deleted.
   
global PGLOBAL;
   
% last argument may be TOL
ni = length(varargin);
if ni>=2,
   lastarg = varargin{ni};
   if isa(lastarg,'char'),
      ni = ni-1; tol = str2double(lastarg);
      if ~isempty(lastarg),
         if isnan(tol) | length(tol)~=1 | ~isreal(tol) | tol<0 | tol>1,
            error('Invalid last argument.');
         end;
      end;     
   end;
end;

% process the 1st input
Pom = [];
eval('Pom = pol(varargin{1});','error(''Invalid argument.'');');
Pdeg = Pom.d;
Pcoef = Pom.c;
Ps = Pom.s;
Pv = Pom.v;
Ph = Pom.h;

% size checking, variable consistency checking, tsp case checking,
if ni>=2,              % sampling periods checking
   for i = 2:ni,
      pom = varargin{i};
      eval('pom = pol(pom);', ...
         'eval(''pom = mdf(pom);'',''error(''''Invalid argument.'''');'');');   
      if isa(pom,'mdf'),             % mdf case
         varargin{i} = pom;
         eval('C = defract(horzcat(varargin{1:ni}));',...
            'error(peel(lasterr));');
         return;
      end;
      if pom.s(1)~=Ps(1),
         error('Matrix has not the same number of rows as the other ones.');
      end;
      if ~isempty(pom.d) & ~isempty(pom.v) & ~strcmp(Pv,'inc'),
         if strcmp(Pv,''),
            Pv=pom.v;
         elseif strcmp(Pv,'z'),
            if strcmp(pom.v,'z^-1'), Pv = 'zz^-1';
            elseif ~strcmp(pom.v,'z'), Pv = 'inc';
            end;
         elseif strcmp(Pv,'z^-1'),
            if strcmp(pom.v,'z'), Pv = 'zz^-1';
            elseif ~strcmp(pom.v,'z^-1'), Pv = 'inc';
            end;
         elseif strcmp(Pv,'zz^-1'),
            if ~(strcmp(pom.v,'z') | strcmp(pom.v,'z^-1')), Pv = 'inc'; nic = i;
            end;
         else,
            if ~(strcmp(pom.v,Pv)), Pv = 'inc';
            end;
         end;
      end;
      if isempty(Ph),
         Ph = pom.h;
      elseif ~isnan(Ph) & Ph~=-1,
         if ~isempty(pom.h),
            if isnan(pom.h),
               Ph = NaN;
            elseif pom.h~=Ph,
               Ph = -1;
            end;
         end;
      end;  
   end;
end;

if strcmp(Pv,'inc'),
   Pv = PGLOBAL.VARIABLE;
   warning('Inconsistent variables.');
elseif strcmp(Pv,'zz^-1'),        % tsp case
   if ~isempty(Ph) & Ph==-1,
      warning('Inconsistent sampling periods.');
   end;
   var = cell(1,ni);
   for i = 1:ni,
      var{i} = tsp(varargin{i});
      if ~isempty(Ph) & (isnan(Ph) | Ph==-1),
         var{i}.h = NaN;
      end;
   end;
   C = horzcat(var{1:ni});
   return;
else
   if ~isempty(Ph) & Ph==-1,
      warning('Inconsistent sampling periods.');
      Ph = NaN;
   end;   
end;
if strcmp(Pv,'s') | strcmp(Pv,'p'),
   Ph = 0;
end;

%concatenation
if ni>=2,
   for i = 2:ni,
      pom = pol(varargin{i});
      if isempty(Pdeg),
         Pdeg = pom.d;
      elseif ~isempty(pom.d),
         Pdeg = max(Pdeg,pom.d);
      end;   
     
      % new coefficient array
      s = size(pom.c,3) - size(Pcoef,3);
      Pcoef = cat( 2,cat(3,Pcoef,zeros([Ps,s])) , ...
                     cat(3,pom.c,zeros([pom.s,-s])) );
         
      % new size
      Ps = [Ps(1), Ps(2) + pom.s(2)]; 
   end;                  
end;

C.d = Pdeg;
C.s = Ps;
C.c = Pcoef;
C.v = Pv;
C.h = Ph;
C.u = [];
C.version = 3.0; 

C = class(C, 'pol');

%end.. @pol/horzcat

