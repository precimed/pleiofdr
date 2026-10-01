function [logpvec1, logpmat2] = correct_sample_overlap(logpvec1, logpmat2, ivec0, traitname1, traitnames)
    % Decorrelation in Statistics: The Mahalanobis Transformation
    % Added material to Data Compression: The Complete Reference
    % http://www.davidsalomon.name/DC2advertis/DeCorr.pdf

    if size(logpmat2, 2) ~= 1, error('unable to control for sample overlap with more than 2 traits\n'); end;

    tmp_zvec1 = -norminv( 0.5 * (10 .^(-logpvec1) ));
    tmp_zmat2 = -norminv( 0.5 * (10 .^(-logpmat2) ));

    idx = ivec0 & ~isnan(tmp_zvec1) & ~isnan(tmp_zmat2);
    C = corrcoef(tmp_zvec1(idx), tmp_zmat2(idx));
    fprintf('Correct correlation between z scores in %s and %s; before correction the correlation was %.3f\n', traitname1, traitnames{1}, C(1,2))

    z1z2 = C^(-0.5) * [tmp_zvec1, tmp_zmat2]';
    tmp_zvec1 = z1z2(1, :)';
    tmp_zmat2 = z1z2(2, :)';

    logpvec1 = -log10(2*normcdf(-abs(tmp_zvec1)));
    logpmat2 = -log10(2*normcdf(-abs(tmp_zmat2)));
end
