function varargout = impulse(F,varargin)
% FRAC/IMPULSE   This function is a Polynomial Toolbox wrapper 
%   for the Control Systems Toolbox function IMPULSE. This 
%   function enables the IMPULSE command to be called with 
%   Polynomial Toolbox objects (fractions) as input 
%   arguments directly, eg. IMPULSE(1/s+1).

if nargout,
    if nargin>1,
        [varargout{1:nargout}] = impulse(ss(F),varargin{:});
    else
        [varargout{1:nargout}] = impulse(ss(F));
    end;
else
    if nargin>1,
        impulse(ss(F),varargin{:});
    else
        impulse(ss(F));
    end;

end;