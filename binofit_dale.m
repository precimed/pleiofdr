% get upper and lower 95% confidence limits by formula from the Engineering

% Statistics Handbook at http://www.itl.nist.gov/div898/handbook/prc/section2/prc241.htm

%

function [phat,cint] = binofit_dale(x,n)

    phat = x./n;

    n2 = n + n;

    z = 1.96; zz = z*z./(n2);

    n1 = phat + zz;

    n2 = z*sqrt((phat.*(1-phat)./n + zz./n2));

    d = 1 + zz + zz;

    % my edit - kg (took out transpose, semicolon)

    % cint = [(n1-n2)./d; (n1+n2)./d]';

    cint = [(n1-n2)./d (n1+n2)./d];

end    

    

