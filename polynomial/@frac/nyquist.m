function varargout = nyquist(F,varargin)
% FRAC/NYQUIST   This function is a Polynomial Toolbox wrapper 
%   for the Control Systems Toolbox function NYQUIST. This 
%   function enables the NYQUIST command to be called with 
%   Polynomial Toolbox objects (fractions) as input 
%   arguments directly, eg. NYQUIST(1/s+1).

if nargout,
    if nargin>1,
        [varargout{1:nargout}] = nyquist(ss(F),varargin{:});
    else
        [varargout{1:nargout}] = nyquist(ss(F));
    end;
else
    if nargin>1,
        nyquist(ss(F),varargin{:});
    else
        nyquist(ss(F));
    end;

end;