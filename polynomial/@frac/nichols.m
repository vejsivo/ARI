function varargout = nichols(F,varargin)
% FRAC/NICHOLS   This function is a Polynomial Toolbox wrapper 
%   for the Control Systems Toolbox function NICHOLS. This 
%   function enables the NICHOLS command to be called with 
%   Polynomial Toolbox objects (fractions) as input 
%   arguments directly, eg. NICHOLS(1/s+1).

if nargout,
    if nargin>1,
        [varargout{1:nargout}] = nichols(ss(F),varargin{:});
    else
        [varargout{1:nargout}] = nichols(ss(F));
    end;
else
    if nargin>1,
        nichols(ss(F),varargin{:});
    else
        nichols(ss(F));
    end;

end;