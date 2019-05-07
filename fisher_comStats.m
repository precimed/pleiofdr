function lp = fisher_comStats(lp1, lp2, ivec0, stdgc, pruneidx, excludevec)
%% FISHER Make Fisher combined stats
   
    if exist('excludevec', 'var')
        lp1(excludevec) = NaN;
        for iteri=1:size(lp2,2), lp2(excludevec, iteri) = NaN; end
    end
    if ~exist( 'pruneidx', 'var' )
        pruneidx = [];
    end
    if ~exist( 'stdgc', 'var' )
        stdgc = false;
    end
    lp = NaN(size(lp1,1), size(lp2,2));
    for iteri = 1:size(lp2,2)
        lp(:, iteri) = -log10(gammainc( ((norminv(10 .^ (-lp1) / 2)) .^ 2 +...
            (norminv(10.^(-lp2(:,iteri))/2)).^2)/2, 2/2, 'upper'));  %????
        if ~exist('ivec0', 'var')
            lp(:,iteri) = GCcorrect_logpvec(lp(:, iteri), ivec0, stdgc, pruneidx);
        end
    end
