function save_zscore(results, traitname1, traitnames, opts, outdir, pruned)
    % 14.08 different names when saving unpruned list
    % 29.08 line 43: corrected mistake A2 = results.A1vec(i);
    
    if ~exist('outdir', 'var'), outdir = '~/Desktop'; end
    if ~exist('pruned', 'var'), pruned=true; end
    ncondtraits = length(traitnames);
    
    if pruned
        fname = sprintf('%s/%s%s_zscore_%s_%g_loci.csv', outdir, traitname1, ...
            sprintf('_%s', traitnames{:}), opts.stattype, opts.fdrthresh);
    else
        fname = sprintf('%s/%s%s_zscore_%s_%g_all.csv', outdir, traitname1, ...
            sprintf('_%s', traitnames{:}), opts.stattype, opts.fdrthresh);
    end
    
    fid = fopen(fname,'w');
    if pruned
        ivec = find(isfinite(results.locusnumvec) & sum(results.imat2, 2) > 0);
    else
        ivec = find(isfinite(results.locusnumvec));
    end
    fprintf(fid, 'locusnum,snpid,geneid,chrnum,chrpos,A1,A2,zscore_%s', traitname1);
    for j = 1:ncondtraits, fprintf(fid, ',zscore_%s', traitnames{j}); end
    for j = 1:ncondtraits 
        fprintf(fid, ',%s_%s_%s,prune_%s_%s', opts.stattype, traitname1, traitnames{j}, ...
            traitname1, traitnames{j});
    end
    fprintf(fid, ',min_%s', opts.stattype);
    fprintf(fid, ',pval_%s', traitname1);
    for j = 1:ncondtraits
        fprintf(fid,',pval_%s', traitnames{j});
    end
    fprintf(fid,'\n');
    for ii = 1:length(ivec) 
        i = ivec(ii);
        z = results.zvec(i);
        A1 = results.A1vec{i};
        A2 = results.A2vec{i};
        fprintf(fid, '%d,%s,%s,%d,%d,%s,%s,%f', results.locusnumvec(i), results.snpidlist{i}, ...
            results.genenamelist{i}, results.chrnumvec(i), results.posvec(i), A1, A2, z);
        for j = 1:ncondtraits, fprintf(fid,',%f', results.zmat2(i, j)); end
        for j = 1:ncondtraits
            fprintf(fid, ',%e,%d', 10.^-results.logfdrmat(i,j), results.imat2(i,j)); 
        end
        fprintf(fid, ',%e', min(10.^-results.logfdrmat(i,:)));
        fprintf(fid, ',%e', 10.^-results.logpvec(i));
        for j = 1:ncondtraits
          fprintf(fid, ',%e', 10.^-results.logpmat2(i,j));
        end
        fprintf(fid, '\n');
    end
    fclose(fid); 
end
