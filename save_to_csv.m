function save_to_csv(results, opts, traitname1, traitnames, outputdir, prunecsv)
% SAVE_TO_CSV Prints results on output file at directory OUTDIR
%    Default prunecsv=true lists only the most significant SNP between linked SNPs
%    Inputs "outdir" and "prunecsv" are optional.

    if ~exist('outputdir', 'var'), outputdir = '~/Desktop';end
    if ~exist('prunecsv', 'var'), prunecsv = true; end
    
    mkdir(outputdir);

    if strcmpi(opts.stattype, 'condfdr')
        save_fdr(results, traitname1, traitnames, opts, outputdir, prunecsv);
    elseif strcmpi(opts.stattype, 'conjfdr')
        save_fdr(results, traitname1, traitnames, opts, outputdir, prunecsv);
        save_zscore(results, traitname1, traitnames, opts, outputdir, prunecsv);
    else
        error ('Not Implemented FDR type %s', opts.stattype);
    end

end
