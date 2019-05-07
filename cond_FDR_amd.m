function [fdrvec, lookup, fdrvec0 ] = cond_FDR_amd(lp1, lp2, opts, pruneidx, excludevec)
% cond_FDR   Recalculate FDR from -log10(p1) condition on -log10(p2)
%      using a lookup table
% 02.07: CRITICAL: [fdrvec, fdrvec0] is now [fdr,lookup, fdrvec0]
%        added lookuptable counts and standard deviation,
%        output as lookup = {mean, count, std}
% 16.01: made pruneidx always an input, change :
%        if exist('pruneidx') to
%        if ~isempty(pruneidx)

% AMD: remove excluded SNPs for purposes of lookup_table creation
lp1_tmp = lp1; lp1_tmp(excludevec,:) = NaN;
lp2_tmp = lp2; lp2_tmp(excludevec,:) = NaN;

if ~isempty(pruneidx)
    if size(pruneidx, 1) ~= size(lp1,1)
        error ('random prune index not conform to logpvec 1');
    end
    if size(pruneidx, 1) ~= size(lp2,1)
        error ('random prune index not conform to logpvec 2');
    end
    [looktable, lookcount, lookstd] = lookup_table(lp1_tmp, lp2_tmp, opts, pruneidx);
else
    [looktable, lookcount, lookstd] = lookup_table(lp1_tmp, lp2_tmp, opts);
end
lookup = {looktable, lookcount, lookstd};
fdrvec0 = NaN;
hv = opts.t1breaks(1:opts.thin:end);

if nargout == 3
    indvec1 = min(length(hv), max(1, 1+(lp1 / (hv(2) - hv(1)))));
    indvec2 = ones(size(indvec1));
    fdrvec0 = interp2(looktable, indvec1, indvec2);
end
indvec1 = min(length(hv),max(1,1 + (lp1 / (hv(2) - hv(1)))));
indvec2 = min(opts.t2_nbreaks, max(1, 1 + (lp2 / (opts.t2breaks(2) - opts.t2breaks(1)))));
fdrvec = interp2(looktable, indvec1, indvec2);
fdrvec(~isfinite(lp1 + lp2)) = NaN;

end
