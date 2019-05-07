function check_sample_overlap(logpvec1, logpmat2, ivec0, traitname1, traitnames)
    % Decorrelation in Statistics: The Mahalanobis Transformation
    % Added material to Data Compression: The Complete Reference
    % http://www.davidsalomon.name/DC2advertis/DeCorr.pdf

    fprintf('Testing for sample overlap...\n');
    tmp_zvec1 = -norminv( 0.5 * (10 .^(-logpvec1) ));
    tmp_zmat2 = -norminv( 0.5 * (10 .^(-logpmat2) ));
    for j=1:size(logpmat2,2)
        idx = ivec0 & ~isnan(tmp_zvec1) & ~isnan(tmp_zmat2(:, j));
        C = corrcoef([tmp_zvec1(idx), tmp_zmat2(idx, j)]);
        fprintf('\tCorrelation between z scores in %s and %s is %.3f\n', ...
            traitname1, traitnames{j}, C(1,2))
    end
    fprintf('\tNote that correlation is calculated across intergenic SNPs. A large correlation\n')
    fprintf('\tmay indicate sample overlap.\n');
end
