function [imat, imat2, logfdrmat] = ind_loci_idx(fdrmat, flp, LDmat, mafvec, opts)
% IND_LOCI_IDX   select SNP indices which passes the p-threshold and
%   FDR-threshold set by options
% 11.06 : changed from flp>-log10(pthresh) to -log10<=pthresh for consistency 

defvec = isfinite(fdrmat + flp);
imat = fdrmat <= opts.fdrthresh & 10.^(-flp) <= opts.pthresh;

logfdrmat = -log10(fdrmat);
logfdrmat(~defvec) = NaN;
if ~isnan(opts.mafthresh), logfdrmat(isnan(mafvec) | mafvec <= opts.mafthresh) = NaN; end;
logfdrmat_pruned = logfdrmat;
logfdrmat_pruned(~imat) = NaN;

imat2 = false(size(imat));
for iteri=1:size(fdrmat,2)
    tmp = FastPrune(logfdrmat_pruned(:, iteri), LDmat);

    % Ensure that only one hit per excluded region survives loci pruning
    for exclude_idx=1:length(opts.exclude_from_fit)
       range = opts.exclude_from_fit{exclude_idx};
       range_indices = false(size(logfdrmat_pruned, 1), 1);
       range_indices(range(1):range(2)) = true;
       range_fdr = logfdrmat_pruned(:, iteri);
       range_fdr(~range_indices) = nan;

       tmp(range_indices) = nan;
       if any(~isnan(range_fdr))
           [~, id] = max(range_fdr);
           tmp(id) = logfdrmat_pruned(id, iteri);
       end
    end

    imat2(:, iteri) = isfinite(tmp);
end

return;
