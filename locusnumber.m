function [locusnumvec, locusnum] = locusnumber(imat, LDmat, mafvec, opts)
% 14.11 critical change locusnumvec(mini:maxi) to locusnumber(idx)
%       where idx = intersect(mini:maxi, ivec)
% 14.11 general comments, change to iivec = any(.)

% choose only snp passing threshold
iivec = any(imat, 2);  % imat = (fdrmat <= opts.fdrthresh)
if ~isnan(opts.mafthresh), iivec(isnan(mafvec) | mafvec <= opts.mafthresh, :) = 0; end;
ivec = find(iivec);

% if nothing passes threshold, exit function
if isempty(ivec) 
    locusnumvec = NaN;
    locusnum = NaN;
    return;
end


%initialize first block represented by locus
locusnum = 0;
locusnumvec = NaN(size(iivec));
mini = ivec(1); maxi = mini;

%repeat until no new blocks are found
while ~isempty(mini)
    
    % add locusnum for each new block
    locusnum = locusnum + 1;
    
    % calculate the extent of the locus (in n-degree LD with mini)    
    while 1
        ii = maxi;
        maxi = find(LDmat(:, ii) & iivec, 1, 'last');
        if maxi <= ii, break; end
    end
    
    % label snps with locusnum, preserve NaNs
    idx = intersect(mini:maxi, ivec);   
    locusnumvec(idx) = locusnum;
    
    % prepare indices for next locus
    mini = min(ivec(ivec > maxi));
    maxi = mini;
end
