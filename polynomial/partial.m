function y = partial(f, flag);
% PARTIAL   Partial fraction expansion.
% 
%  The command
%     E = PARTIAL(F)
%  computes the vector E containing the partial fractions 
%  of a scalar real rational function F. F = SUM(E) holds. 
%  E is of the class MDF regardless of the class of F.
%
%  Real or complex decomposition can be selected 
%  by an additional input argument:
%     E = PARTIAL(F,'real')     (default)
%     E = PARTIAL(F,'complex')

% M. Hromcik, August 29, 2005
% September 6


if any(f.s>1), 
    error('Partial fraction expansion works for real scalar fractions only.'); 
end;

if ~isreal(f), 
    error('Partial fraction expansion works for real scalar fractions only.'); 
end;

if nargin == 1, flag='real'; end

var = f.var;
[R,P,K] = residue(pol2mat(f.num), pol2mat(f.den));
y = mdf(zeros(1,length(P)));

switch flag(1),

    case 'c',
    % complex, degree 1 
    
    % real poles first
    indreal = find(~imag(P));
    Preal = P(indreal);
    Rreal = R(indreal); 
    
    indimag = find(imag(P));
    Pimag = P(indimag);
    Rimag = R(indimag);
    
    R = [Rreal;Rimag];
    P = [Preal;Pimag];
    
    i = 1;
    while i<=length(P),    
         j = 0;
         while (i+j)<=length(P) & (P(i) == P(i+j)),
             y(i+j) = mdf(R(i+j),(eval(var)-P(i+j))^(j+1));
             j = j+1;
         end;
         i = i+j;
    end; 
    
    case 'r',
    % real, degree 2

    % separate real and imaginary poles with positive imaginary part
    indreal = find(~imag(P));
    Preal = P(indreal);
    Rreal = R(indreal); 
    
    indimag = find(imag(P)>0);
    Pimag = P(indimag);
    Rimag = R(indimag);

    % process real poles
    i = 1;
    while i<=length(Preal),    
         j = 0;
         while (i+j)<=length(Preal) & (Preal(i) == Preal(i+j)),
             y(i+j) = mdf(Rreal(i+j),(eval(var)-Preal(i+j))^(j+1));
              j = j+1;
         end;
         i = i+j;
    end;

    % process imaginary poles
    shft = i-1;
    i = 1;
    while i<=length(Pimag),    
         j = 0;
         while (i+j)<=length(Pimag) & (Pimag(i) == Pimag(i+j)),
             %y(shft+i+j) = mdf( 2*real(Rimag(i+j))*s - 2*real(Rimag(i+j))*real(Pimag(i+j))-2*imag(Rimag(i+j))*imag(Pimag(i+j)),(eval(var)^2 - 2*real(Pimag(i+j))*eval(var) + abs(Pimag(i+j))^2)^(j+1) );
             y(shft+i+j) = real( mdf(Rimag(i+j),(eval(var)-Pimag(i+j))^(j+1)) + mdf( conj(Rimag(i+j)),(eval(var)-conj(Pimag(i+j)))^(j+1)) );
             j = j+1;
         end;
         i = i+j;
    end;
    
    if any(y==0),
        y = y(:,1:min(find(y==0))-1);
    end;

end; % switch

if ~isempty(K), 
    polpart = lop(K, length(K)-1); 
    polpart.v = var;
    y = [mdf(polpart,1) y];
end;



