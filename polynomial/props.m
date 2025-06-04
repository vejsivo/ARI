function props(P,varargin)
%PROPS  Display/modify properties of object
%
% The commmand
%   PROPS(P,VALUE)  
% sets the property of P corresponding to VALUE equal to the value 
% VALUE. The commmand
%   PROPS(P,Value1,Value2,...)  
% sets multiple property values of P equal to the values Value1, 
% Value2,... The commmand
%   PROPS(P)  
% displays all the properties of P and their admissible values.
%
% See also: POL/PROPS, TSP/PROPS, FRAC/PROPS.

%       Author(s):  S. Pejchova  16-03-2000
%       Copyright (c) 1999 by Polyx, Ltd.
%       $Revision: 3.0 $  $Date:16-Mar-2000 11:12:50   $
%                         $Date:28-Feb-2003 $

ni = nargin;

if ni==1,
  disp(' ');
  disp('The input is not a polynomial object.');
  disp(['It is a ',cl2str(P),'.']);
  disp(' ');
else,
   error('Use SET to modify the properties of the nonpolynomial object.');
end

%end .. props
