function [logpvec_gc, sig0] = GCcorrect_logpvec(logpvec,ivec0,use_standard_gc,pruneidx)
% GCCORRECT_LOGPVEC  Perform a correction for genomic control on -log(p)
%
% Reference: Devlin, B., and Kathryn Roeder. "Genomic control for
%    association studies." Biometrics 55.4 (1999): 997-1004.

if ~islogical(ivec0), error('ivec0 must be a logical vector of control elements'); end
if ~exist('use_standard_gc', 'var'), use_standard_gc = false; end
if ~exist('pruneidx', 'var'), pruneidx = []; end

if isempty(pruneidx), pruneidx = true(size(logpvec)); end

randprune_n = size(pruneidx, 2);
sig0 = nan(1, randprune_n);

for randprune_i = 1:randprune_n

    % Get median
    logpvec0 = logpvec(ivec0 & pruneidx(:, randprune_i));
    zvec0 = abs(norminv(1/2*10 .^ -logpvec0));

    if use_standard_gc
        sig0(randprune_i) = sqrt(median(zvec0(~isnan(zvec0)).^2) ./ chi2inv(0.5,1));
    else
        fracvec = 1 - logspace(log10(0.5), -2, 50);
        prc = prctile(zvec0.^2, 100*fracvec);
        prc = reshape(prc, 1, length(prc)); % octave-fix: convert prc nx1 to 1xn
        sig0vec = sqrt(prc ./ chi2inv(fracvec, 1));
        sig0(randprune_i) = median(sig0vec);
    end

end

% Normalize by median
zvec = abs(norminv(1/2*10 .^ -logpvec));
logpvec_gc = -log10(2*normcdf(-abs(zvec/median(sig0))));

end

