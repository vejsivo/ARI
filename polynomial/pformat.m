function pformat(string),
%PFORMAT  Set the output format for polynomial matrices
%
% Used to switch between different output display formats 
% as follows:
%    PFORMAT         Default. Same as PFORMAT NICE
%    PFORMAT NICE    Pretty and very simple symbolic form of the polynomial matrix
%    PFORMAT COEF    Coefficient matrices according to increasing powers
%    PFORMAT RCOEF   Coefficient matrices according to decreasing powers
%    PFORMAT BLOCK   Row block coefficient matrix
%    PFORMAT SYMB    Symbolic form of the polynomial matrix
%    PFORMAT SYMBS   Short symbolic form of the polynomial matrix
%    PFORMAT SYMBR   Rational symbolic form of the polynomial matrix
%    PFORMAT ROOTC   Zero/pole/gain form of the polynomial matrix
%    PFORMAT ROOTR   Real zero/pole/gain form of the polynomial matrix
%    PFORMAT NORMAL  According to increasing powers
%    PFORMAT REVERSE According to decreasing powers

%       Author(s):  S. Pejchova, M. Sebek 3-3-98
%       Copyright (c) 1998 by Polyx, Ltd.
%       $Revision: 2.0 $  $Date: 24-Sep-1998 15:18:34   $
%       $Revision: 3.0 $  $Date: 09-May-2000 16:29:00  S.Pejchova  $
%                         $Date: 06-Sep-2000  S.Pejchova  $
%                         $Date: 23-Sep-2001  J.Jezek     $
%                         $Date: 28-Feb-2003  J.Jezek     $

global PGLOBAL;
eval('PGLOBAL.FORMAT;', 'painit;');

if nargin == 0,
   PGLOBAL.FORMAT = 'nice';
   PGLOBAL.ORDER = 'normal';
elseif ~ischar(string), 
   error('Invalid argument; must be a string.');
else,
   switch lower(string),
   case {'symbs','symbr','symb','coef','rcoef','block','rootr','rootc','nice'},
      PGLOBAL.FORMAT = lower(string);
   case {'normal','reverse'},
      PGLOBAL.ORDER = lower(string);
   otherwise,
      error('Invalid format.');
   end;
end;

%end .. pformat
