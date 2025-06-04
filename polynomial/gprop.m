function gprop(varargin)
%GPROP  Set global polynomial propertis
%
% GPROP VALUE sets the value of the related global property
% to VALUE. Works with any number of input arguments such as
% GPROP VALUE1 VALUE2 ... VALUEN. When called without arguments, 
% GPROP displays the current values of all the global properties 
% along with the list of alternatives. 
  
%       Author(s):  S. Pejchova, M. Sebek 3-3-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 10-Nov-1998 14:17:50   $
%       $Revision: 3.0 $  $Date: 16-Feb-2000 12:15:00  J.Jezek  $ 

ni = nargin;
gprops(varargin{1:ni});
%for compatibility, GPROP is the same as GPROPS
        
%end .. gprop
