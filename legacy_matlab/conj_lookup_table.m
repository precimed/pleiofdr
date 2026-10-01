function looktable = conj_lookup_table(lp1, lp2, opts, pruneidx);
    if exist('pruneidx', 'var')
        look12 = lookup_table(lp1, lp2, opts, pruneidx);
        look21 = lookup_table(lp1, lp2, opts, pruneidx);
    else
        look12 = lookup_table(lp1, lp2, opts);
        look21 = lookup_table(lp1, lp2, opts);
    end
    
    [M,N] = size(look12);
    look12 = [ look12 ; repmat(look12(end,:), [N-M, 1]) ];
    look21 = [ look21 ; repmat(look21(end,:), [N-M, 1]) ];
    looktable = max( look12,look21' );
    
%    looktable = max( ...
%        cat(1, look12, repmat(look12(end,:),[size(look12,2)-size(look12,1) 1]) ),...
%        cat(1, look21, repmat(look21(end,:),[size(look21,2)-size(look21,1) 1]) )');
