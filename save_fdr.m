function save_fdr(results, traitname1, traitnames, opts, outdir, pruned)
% 02.07 fixed bug for length of traitnames (line 21)
% 14.08 different names when saving unpruned list
    
    ncondtraits = length(traitnames);
    
    if ~exist('outdir', 'var'), outdir = '/tmp'; end
    if ~exist('pruned', 'var'), pruned=true; end
    
    if pruned
        fname = sprintf('%s/%s%s_%s_%g_loci.csv', outdir, traitname1, ...
            sprintf('_%s', traitnames{:}), opts.stattype, opts.fdrthresh);
    else
        fname = sprintf('%s/%s%s_%s_%g_all.csv', outdir, traitname1, ...
            sprintf('_%s', traitnames{:}), opts.stattype, opts.fdrthresh);
    end
    fid = fopen(fname,'w'); 
    if pruned
        ivec = find(isfinite(results.locusnumvec) & sum(results.imat2,2)>0);
    else
        ivec = find(isfinite(results.locusnumvec));   %AW: ?
    end
    
    %Header
    fprintf(fid, 'locusnum,snpid,geneid,chrnum,chrpos,pval_%s,fdr_%s', traitname1, traitname1);
    for j = 1:ncondtraits
        %fprintf(fid, ',pcomb_%s_%s,%s_%s_%s,prune_%s_%s', traitname1, traitnames{j}, ...
        %    opts.stattype, traitname1, traitnames{j}, traitname1, traitnames{j}); 
        fprintf(fid, ',%s_%s_%s,prune_%s_%s', opts.stattype, traitname1, traitnames{j}, ...
            traitname1, traitnames{j});
    end
    fprintf(fid, ',min_%s', opts.stattype); fprintf(fid,'\n');
    
    %Entries
    for ii = 1:length(ivec)
        i = ivec(ii);
        fprintf(fid, '%d,%s,%s,%d,%d,%e', results.locusnumvec(i), results.snpidlist{i}, ...
            results.genenamelist{i}, results.chrnumvec(i), results.posvec(i), ...
            10 .^ -results.logpvec(i)); 
        fprintf(fid, ',%e', results.fdrvec0(i)); 
        for j = 1:ncondtraits 
            %fprintf(fid, ',%e,%e,%d', 10.^-(results.flp(i,j)), results.fdrmat(i,j), ...
            %    results.imat2(i,j));
            fprintf(fid, ',%e,%d', results.fdrmat(i,j), ...
                results.imat2(i,j)); 
        end 
        fprintf(fid, ',%e', min(10 .^ -results.logfdrmat(i,:))); 
        fprintf(fid, '\n'); 
    end 
    fclose(fid); 
end
