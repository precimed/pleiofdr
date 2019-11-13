% (c) 2019 University of Oslo
% Pleiotropy-informed conditional and conjunctional false discovery rate

% Running
%   - pleioFDR         : Conditional / Conjunction FDR
%   - plotqq           : Q-Q plots
%   - plot_enrichment  : Fold Enrichment
%   - plot_manhattan   : Manhattan plots
%
% Requires the following data to be in your workspace:
%
%   options        : an instance of pleioOpt class with options and settings for this run
%   traitfile1     : string, .mat file containing main phenotype
%   traitfiles     : cell array of strings, .mat files for secondary phenotypes
%   traitname1     : string, name of the main phenotype
%   traitnames     : cell array of strings, name of secondary phenotypes
%   LDmat          : logical sparse matrix, nsnp x nsnp, LD matrix cut at certain threshold
%   chrnumvec      : matrix, nsnp x 1, chromosome (1..22)
%   posvec         : matrix, nsnp x 1, base pair position on chromosome
%   ivec0          : logical matrix, nsnp x 1, indicating variants to use for genomic correction
%                    (e.g. intergenic)

%% LOAD OCTAVE PACKAGES

if is_octave()
    pkg load statistics;
    pkg load nan;
    graphics_toolkit('gnuplot');
end


%% INPUT VALIDATION

if ~exist('options', 'var'), options = pleioOpt(); end
if ~isempty(options.outputdir), outputdir = options.outputdir; end
if ~exist('outputdir', 'var') || isempty(outputdir), outputdir='test'; end
if exist('fishercomb', 'var') && fishercomb ~= options.fishercomb
        error('fishercomb must be specified via pleioOpt options');
end
if exist('correct_for_sample_overlap', 'var') && ...
        correct_for_sample_overlap ~= options.correct_for_sample_overlap
    error('correct_for_sample_overlap must be specified via pleioOpt options');
end


%% START DIARY

mkdir(outputdir)
diary(sprintf('%s/%s%s_%s_%g.log', outputdir,traitname1,...
        sprintf('_%s',traitnames{:}), options.stattype, options.fdrthresh))
display(options)



%% LOAD FILES

fprintf('Loading GWAS .mat files... ')
[logpvec1,zvec1]=load_gwas(traitfile1, length(chrnumvec));
[logpmat2,zmat2]=load_gwas(traitfiles, length(chrnumvec));
fprintf('done\n')

%% EXCLUDE SPECIAL SNPS FROM ANALYSIS

excludevec = false(size(logpvec1,1),1);
fprintf('%i variants are defined across all traits\n', sum(isfinite(logpvec1+sum(logpmat2,2))));
if ~isnan(options.mafthresh)
    fprintf('Exclude %i variants due to MAF (minor allele frequency) below %.3f or undefined\n', ...
        sum(isnan(mafvec) | (mafvec <= options.mafthresh)), options.mafthresh)
	excludevec = excludevec | isnan(mafvec) | (mafvec <= options.mafthresh);
end

for j = 1:length(options.exclude_from_fit)
    
     % Should only exclude for fdr calculations, not discovery
    excludevec(options.exclude_from_fit{j}(1):options.exclude_from_fit{j}(2)) = true;

    exclude_from = options.exclude_from_fit{j}(1);
    exclude_to = options.exclude_from_fit{j}(2);
    exclude_chr = unique(chrnumvec(exclude_from:exclude_to));

    fprintf('Exclude %i SNPs on chromosome %s from %i to %i KB\n', ...
        exclude_to - exclude_from + 1, ...
        mat2str(exclude_chr), posvec(exclude_from)/1000, posvec(exclude_to)/1000)
end

excludevec_discovery = false(size(excludevec));
if options.exclude_from_fit_and_discovery
    fprintf('Excluded SNPs will be excluded from both fit and discovery')
    excludevec_discovery = excludevec;
end

if options.exclude_ambiguous_snps
    excludevec = excludevec | is_ambiguous;
    excludevec_discovery = excludevec_discovery | is_ambiguous;
    fprintf('%i Ambigous SNPs will be excluded from both fit and discovery\n', sum(is_ambiguous))
end

if any(excludevec_discovery)
    logpmat2(excludevec_discovery, :) = NaN;
    logpvec1(excludevec_discovery) = NaN;

    zmat2(excludevec_discovery, :) = NaN;
    zvec1(excludevec_discovery) = NaN;
end

if ~any(excludevec), fprintf('No variants were excluded'); end
if any(excludevec), fprintf('%i variants left after exclusion described above\n', ...
        sum(~excludevec & isfinite(logpvec1+sum(logpmat2,2))));
end

%% RANDOM PRUNING INDICES

if options.reset_pruneidx || ~exist('pruneidx', 'var'), pruneidx = []; end
if options.randprune && isempty(pruneidx)
    % AMD: make sure we only consider variants with defined values when pruning
    defvec = ~excludevec & isfinite(logpvec1+sum(logpmat2,2));
    pruneidx = random_prune_idx_amd_fb(options.randprune_n, LDmat, defvec, options.randprune_repeats);
end
if options.randprune && (size(pruneidx, 2) ~= options.randprune_n)
    error('Invalid pruneidx; try reset_pruneidx=false.');
end

%% GENOMIC CONTROL

if(options.perform_gc)
    fprintf('\nPerforming genomic correction... \n')
    fprintf('Use %i control variants to calculate lambda GC (genomic correction factor)\n', ...
        sum(ivec0 & ~excludevec));
    if options.randprune_gc && options.randprune
        pruneidx_gc = pruneidx;
        fprintf('Random pruning is taken into account in lambdaGC calculation. \n')
    else
        pruneidx_gc = [];
    end
    [logpvec1, sig0] = GCcorrect_logpvec(logpvec1, ivec0, options.use_standard_gc, pruneidx_gc);
    fprintf('Genomic Control for %s: lambda = %.3f\n', traitname1, median(sig0));
    zvec1 = -norminv(10.^-logpvec1/2).*sign(zvec1);
    for i=1:length(traitfiles)
        [logpmat2(:,i), sig0] = GCcorrect_logpvec(logpmat2(:,i), ivec0, options.use_standard_gc, ...
            pruneidx_gc);
        fprintf('Genomic Control for %s: lambda = %.3f\n', traitnames{i}, median(sig0));
        zmat2(:,i)    = -norminv(10.^-logpmat2(:,i)/2).*sign(zmat2(:,i));
    end
    fprintf('done\n')
else
    fprintf('Skipping genomic correction. \n')
end


%% CORRECT FOR SAMPLE OVERLAP
% Decorrelation in Statistics: The Mahalanobis Transformation
% Added material to Data Compression: The Complete Reference
% http://www.davidsalomon.name/DC2advertis/DeCorr.pdf

% Calculation of 'FDRvec0' and 'save_to_csv' should use the original (non-decorrelated) stats
% why? and does this comment actually belong here? -NoF
logpvec1_orig = logpvec1; zvec1_orig = zvec1;
logpmat2_orig = logpmat2; zmat2_orig = zmat2;
if options.correct_for_sample_overlap
    [logpvec1, logpmat2] = correct_sample_overlap(logpvec1, logpmat2, ivec0, traitname1, traitnames);
    zvec1 = -norminv( 0.5 * (10 .^(-logpvec1) )) .* sign(zvec1);
    zmat2 = -norminv( 0.5 * (10 .^(-logpmat2) )) .* sign(zmat2);
else
    check_sample_overlap(logpvec1, logpmat2, ivec0, traitname1, traitnames);
end

%% FISHER COMBINED STATISTICS

fprintf('Calculating Fisher combined stats... ')
if (options.fishercomb)
    if options.randprune_gc && options.randprune
        pruneidx_fc = pruneidx;
    elseif options.randprune
        pruneidx_fc = [];
    end
% AMD: Should not exclude exclude_from_fit variants here; should distinguish from exclude_from_plot!
%     flp = fisher_comStats(logpvec1, logpmat2, ivec0, options.use_standard_gc, pruneidx_fc, excludevec);
    flp = fisher_comStats(logpvec1, logpmat2, ivec0, options.use_standard_gc, pruneidx_fc);
else
    flp = repmat(logpvec1, [1, size(logpmat2,2)]);
end
for iteri=1:size(logpmat2,2)
    exclude = isnan(logpmat2(:,iteri));
    flp( exclude, iteri) = NaN;
end

fprintf('done\n')

%% PLEIOTROPY ANALYSIS

fprintf('Running pleiotropy analyses... ')
[fdrmat, ~, ~, fdrmat12, fdrmat21, lookup12, lookup21] = ...
    pleioFDR_amd(logpvec1, logpmat2, options, LDmat, excludevec, pruneidx);
% the unconditioned fdr is calculated from unpruned variants
[~ , ~, fdrvec0] = cond_FDR_amd(logpvec1_orig, logpmat2(:,1), options, [], []);
fprintf('done\n')
fprintf('mean (std) variants per random pruning iteration = %.2f (%.2f)\n', mean(sum(pruneidx)), ...
    std(sum(pruneidx)))

%% PLOT QQ AND ENRICHMENT

fprintf(1,'PLOT QQ / Enrichment... ')
[h_qq, mat_qq]                 = plot_qq_amd(logpvec1, logpmat2, traitname1, traitnames, options, LDmat, pruneidx, excludevec);
[h_qq_inv, mat_qq_inv]         = plot_qq_amd(logpvec1, logpmat2, traitname1, traitnames, options, LDmat, pruneidx, excludevec, true);
[h_tdr, mat_tdr]               = plot_qq_amd(logpvec1, logpmat2, traitname1, traitnames, options, LDmat, pruneidx, excludevec, false, false, true);
[h_tdr_inv, mat_tdr_inv]       = plot_qq_amd(logpvec1, logpmat2, traitname1, traitnames, options, LDmat, pruneidx, excludevec, true, false, true);
[h_enrich, mat_enrich]         = plot_enrichment_amd(logpvec1,logpmat2, traitname1, traitnames, options, LDmat, pruneidx, excludevec);
[h_enrich_inv, mat_enrich_inv] = plot_enrichment_amd(logpvec1,logpmat2, traitname1, traitnames, options, LDmat, pruneidx, excludevec, true);
qq_t1breaks                    = options.t1breaks;

save(fullfile(options.outputdir, 'result.mat'), '-v6', 'fdrmat', 'logpvec1', 'logpmat2', 'zvec1', 'zmat2', 'excludevec', 'mafvec', 'mat_qq', 'mat_qq_inv', 'mat_enrich', 'mat_enrich_inv', 'qq_t1breaks');

mkdir(outputdir)
filetypes = {'fig','jpg','png'};
for j = 1:length(filetypes)
    save_figure(h_qq,         filetypes{j}, sprintf('%s/%s_vs_%sqq', outputdir, ...
        traitname1, sprintf('%s_',traitnames{:})));
    save_figure(h_qq_inv,     filetypes{j}, sprintf('%s/%svs_%s_qq', outputdir, ...
        sprintf('%s_',traitnames{:}), traitname1));
    save_figure(h_tdr,        filetypes{j}, sprintf('%s/%s_vs_%stdr', outputdir, ...
        traitname1, sprintf('%s_',traitnames{:})));
    save_figure(h_tdr_inv,    filetypes{j}, sprintf('%s/%svs_%s_tdr', outputdir, ...
        sprintf('%s_',traitnames{:}), traitname1));
    save_figure(h_enrich,     filetypes{j}, sprintf('%s/%s_vs_%senrich', outputdir, ...
        traitname1, sprintf('%s_',traitnames{:})));
    save_figure(h_enrich_inv, filetypes{j}, sprintf('%s/%svs_%s_enrich', outputdir, ...
        sprintf('%s_',traitnames{:}), traitname1));
end
fprintf('done\n')

%% PLOT LOOKUP

% AMD: how does plot_lookup use logpvec1 & logpmat2 vs. lookup12 and lookup21?
h_lookup    = plot_lookup(logpvec1, logpmat2, traitname1, traitnames, options, LDmat, pruneidx, ...
    excludevec, lookup12, lookup21);

filetypes = {'fig','jpg','png'};
for j = 1:length(filetypes)
    save_figure(h_lookup, filetypes{j}, sprintf('%s/%svs_%s_lookup', outputdir, ...
        sprintf('%s_',traitnames{:}), traitname1));
end

%% SAVE FDR-TABLES AND PLOT MANHATTAN AT A SPECIFIC THRESHOLD

if ~exist('snpidlist', 'var'), snpidlist = cell(size(chrnumvec)); snpidlist(:) = {''}; end;
genenamelist = cell(size(chrnumvec)); genenamelist(:) = {''};
if ~exist('A1vec', 'var'), A1vec = cell(size(chrnumvec)); A1vec(:) = {''}; end;
if ~exist('A2vec', 'var'), A2vec = cell(size(chrnumvec)); A2vec(:) = {''}; end;

% FIND LOCI AT GIVEN THRESHOLD
[imat, imat2, logfdrmat] = ind_loci_idx(fdrmat, flp, LDmat, mafvec, options);
[locusnumvec, locusnum] = locusnumber(imat, LDmat, mafvec, options);
results = struct('fdrmat', fdrmat, 'imat', imat, 'imat2', imat2, ...
    'logfdrmat', logfdrmat, 'locusnumvec', locusnumvec, 'logpmat2', logpmat2_orig, ...
    'logpvec', logpvec1_orig, 'zvec', zvec1_orig, 'zmat2', zmat2_orig, 'fdrvec0', fdrvec0, ...
    'snpidlist', {snpidlist}, 'genenamelist', {genenamelist}, ...
    'chrnumvec', chrnumvec, 'posvec', posvec, 'A1vec', {A1vec}, 'A2vec', {A2vec}, 'flp', flp, ...
    'pruneidx', pruneidx);
if (~is_octave())
    disp(results)
end

% SAVE TABLES
fprintf('Saving .csv... ')
save_to_csv(results, options, traitname1, traitnames, outputdir, false); % don't prune CSV
save_to_csv(results, options, traitname1, traitnames, outputdir, true);  % prune CSV
fprintf('done\n')

% PLOT MANHATTAN
if options.manh_plot
    fprintf('Creating Manhattan plots... ')
    h_Manhattan = plot_Manhattan(results, traitname1, traitnames, chrnumvec, options);

    filetypes = {'fig'};
    for j = 1:length(filetypes)
        save_figure(h_Manhattan, filetypes{j}, sprintf('%s/%s%s_%s_%g_manhattan', outputdir, ...
            traitname1, sprintf('_%s',traitnames{:}), options.stattype, options.fdrthresh));
    end
    fprintf('done\n')
end

%% CLOSE DIARY

diary off
