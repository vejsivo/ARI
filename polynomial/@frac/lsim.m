function varargout = lsim(F,varargin)
% FRAC/LSIM   This function is a Polynomial Toolbox wrapper 
%   for the Control Systems Toolbox function LSIM. This 
%   function enables the LSIM command to be called with 
%   Polynomial Toolbox objects (fractions) as input 
%   arguments directly, eg. LSIM(1/s+1).

if nargout,
    if nargin>1,
        [varargout{1:nargout}] = lsim(ss(F),varargin{:});
    else
        [varargout{1:nargout}] = lsim(ss(F));
    end;
else
    if nargin>1,
        lsim(ss(F),varargin{:});
    else
        lsim(ss(F));
    end;
end;