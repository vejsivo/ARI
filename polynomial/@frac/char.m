function [Cn,Cd] = char(F,n)
%CHAR  Convert a polynomial fraction to strings or to cell arrays of strings.
%
% The commmand
%    [Cn,Cd] = CHAR(F) 
% converts the fraction F containing the scalar polynomial numerator into 
% a string Cn and scalar plynomial denominator into a string Cd. If F.n 
% (or F.d) is a polynomial matrix then the result is a cell array of the same 
% size consisting of strings that correspond to the polynomial elements of F.n
% (or F.d). Each scalar coefficient is taken with respect to the basic Matlab 
% format. The commmand
%    [Cn,Cd] = CHAR(F,N) 
% converts each scalar coefficient into a string representation with 
% a maximum N digits of precision. The commmand
%    [Cn,Cd] = CHAR(F,FORMAT) 
% works like CHAR(F,N) but uses the format string FORMAT for each scalar
% coefficient (see SPRINTF for details). The commmand
%    [Cn,Cd] = CHAR(F,'RAT') 
% is a special case of CHAR(F,FORMAT) and uses rational aproximation of 
% the coefficients.
%    [Cn,Cd] = CHAR(F,'ROOTC') 
% uses zero/pole/gain representation of polynomial entries 
%    [Cn,Cd] = CHAR(F,'ROOTR') 
% is a pretty-form of CHAR(F,'ROOTC') without complex numbers (the complex 
% pairs are multiplied).


%      Author(s): S. Pejchova  07-8-00
%      Copyright (c) 2000 by Polyx, Ltd.
%      $Revision: 3.0 $  $Date: 07-Aug-2000 S. Pejchova   $
%                        $Date: 14-Oct-2002 J. Jezek      $

% Effect on other properties:
% Cn and Cd are standard Matlab cell arrays of strings.

ni = nargin;
if ~ni,
   error('Not enough input arguments.');
elseif ~isa(F,'frac'),
   error('1st argument is not polynomial fraction.');
end;
if ni==1,
   Cn = char(F.num);    Cd = char(F.den);
else,
   Cn ='';
   eval('Cn = char(F.num,n);','error(peel(lasterr));');    
   Cd = char(F.den,n);
end;

%end .. @frac/char
