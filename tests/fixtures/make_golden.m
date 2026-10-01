% (c) 2019 University of Oslo
% Generate golden fixtures for the Python port of pleioFDR.
%
% Run from the repository root with the demo data unpacked there
% (pleioFDR_demo_data.tar.gz):
%
%   matlab -batch "run('tests/fixtures/make_golden.m')"
%
% Writes, under tests/fixtures/:
%   golden_<run>.mat     selected workspace variables after runme.m (inputs, prune
%                        indices, every intermediate array and final result)
%   matlab_<run>/*.csv   the CSV tables written by save_to_csv
%   matlab_<run>/*.png   figures, for visual comparison only
%   unit_matlab.mat      inputs/outputs of individual functions on synthetic data

% run() changes into this script's folder; work from the repository root instead
golden_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
cd(golden_root);
addpath(fullfile(golden_root, 'legacy_matlab'));
golden_dir = fullfile('tests', 'fixtures');

golden_base = {
    'reffile=ref_1kgPhase3eur_LDr2p1_DEMO.mat'
    'traitfolder=.'
    'traitfile1=CTG_COG_2018_DEMO.mat'
    'traitname1=COGchr21'
    'traitfiles={''SSGAC_EDU_2016_DEMO.mat''}'
    'traitnames={''EDUchr21''}'
    'randprune=true'
    'randprune_n=20'
    'exclude_chr_pos=[]'
    'manh_colorlist=[1 0 0; 1 0.5 0 ; 0 0.75 0.75; 0 0.5 0; 0.75 0 0.75; 0 0 1; 0 1 0; 0 1 1]'
    'reset_pruneidx=true'
    'randprune_repeats=default'
    'pthresh=1'
    'perform_gc=true'
    'use_standard_gc=false'
    'randprune_gc=true'
    'exclude_from_discovery=false'
    'mafthresh=0.005'
    'exclude_ambiguous_snps=true'
    'onscreen=false'
    'exit_matlab_upon_completion=false'
};

% name, overrides (later lines win in TextConfig)
golden_runs = {
    'conjfdr', {'stattype=conjfdr', 'fdrthresh=0.05'}
    'condfdr', {'stattype=condfdr', 'fdrthresh=0.01'}
    'condfdr_excl', {'stattype=condfdr', 'fdrthresh=0.05', ...
        'exclude_chr_pos=[21 30000000 32000000]', 'exclude_from_discovery=true', ...
        'use_standard_gc=true', 'randprune_gc=false', 'exclude_ambiguous_snps=false', ...
        'randprune_repeats=maxout'}
};

golden_vars = {'pruneidx', 'logpvec1', 'logpmat2', 'zvec1', 'zmat2', 'excludevec', ...
    'flp', 'fdrmat', 'fdrmat12', 'fdrmat21', 'lookup12', 'lookup21', 'fdrvec0', ...
    'mat_qq', 'mat_qq_inv', 'mat_tdr', 'mat_tdr_inv', 'mat_enrich', 'mat_enrich_inv', ...
    'imat', 'imat2', 'logfdrmat', 'locusnumvec', 'locusnum', 'exclude_chr_pos'};

for golden_k = 1:size(golden_runs, 1)
    clearvars -except golden_* LDmat chrnumvec posvec mafvec is_ambiguous is_intergenic ivec0
    close all
    golden_name = golden_runs{golden_k, 1};
    golden_out = fullfile(golden_dir, ['matlab_' golden_name]);
    config = fullfile(tempdir, ['pleiofdr_golden_' golden_name '.txt']);
    golden_fid = fopen(config, 'w');
    fprintf(golden_fid, '%s\n', golden_base{:}, golden_runs{golden_k, 2}{:}, ...
        ['outputdir=' golden_out]);
    fclose(golden_fid);

    rng(20190507 + golden_k);
    runme

    golden_present = golden_vars(ismember(golden_vars, who));
    save(fullfile(golden_dir, ['golden_' golden_name '.mat']), '-v7', golden_present{:});
    delete(fullfile(golden_out, '*.fig'));
    delete(fullfile(golden_out, '*.svg'));
    delete(fullfile(golden_out, '*.mat'));
end

%% Unit fixtures on synthetic data
clearvars -except golden_dir LDmat
close all
rng(42);
u = struct();

% binofit (Clopper-Pearson, MATLAB Statistics Toolbox)
u.binofit_x = [1; 5; 10; 50; 100; 3];
u.binofit_n = [10; 10; 100; 100; 100; 3];
[u.binofit_phat, u.binofit_pci] = binofit(u.binofit_x, u.binofit_n);

% histc / hist3 conventions as used by lookup_table and plot_qq_amd
u.hist_x = [0; 0.005; 0.01; 0.0149; 0.015; 2.5; 29.99; 30; 31; NaN; Inf];
u.histc_counts = histc(u.hist_x, linspace(0, 30, 3001));

% SparseSmooth2d
u.smooth_im = randn(31, 301);
u.smooth_wim = rand(31, 301);
u.smooth_wim(:, 250:end) = 0;
u.smooth_out = SparseSmooth2d(u.smooth_im, u.smooth_wim, [1e2 1e2]);

% Synthetic -log10(p) with enrichment, some NaN
nsnp = size(LDmat, 1);
lp2 = -log10(rand(nsnp, 1));
lp1 = -log10(rand(nsnp, 1) .^ (1 + 0.5 * (lp2 > 2)));
lp1(rand(nsnp, 1) < 0.01) = NaN;
lp2(rand(nsnp, 1) < 0.01) = NaN;
lp1(1:5) = [0; 1e-12; 29.9; 35; Inf];
u.lp1 = lp1; u.lp2 = lp2;

opts = pleioOpt();
opts.randprune_n = 5;
defvec = isfinite(lp1 + lp2);
u.defvec = defvec;
u.fastprune_in = -log10(rand(nsnp, 1)); u.fastprune_in(~defvec) = NaN;
u.fastprune_out = FastPrune(u.fastprune_in, LDmat);
u.pruneidx = random_prune_idx_amd_fb(opts.randprune_n, LDmat, defvec, 'default');

[u.look_nopr, ~, ~] = lookup_table(lp1, lp2, opts);
[u.look_pr, u.lookcount_pr, u.lookstd_pr] = lookup_table(lp1, lp2, opts, u.pruneidx);
[u.cond_fdrvec, ~, u.cond_fdrvec0] = cond_FDR_amd(lp1, lp2, opts, u.pruneidx, false(nsnp, 1));

ivec0 = rand(nsnp, 1) < 0.3;
u.ivec0 = ivec0;
[u.gc_inhouse, u.gc_inhouse_sig0] = GCcorrect_logpvec(lp1, ivec0, false, []);
[u.gc_std, u.gc_std_sig0] = GCcorrect_logpvec(lp1, ivec0, true, []);
[u.gc_pruned, u.gc_pruned_sig0] = GCcorrect_logpvec(lp1, ivec0, false, u.pruneidx);
u.fisher = fisher_comStats(lp1, lp2, ivec0, false, []);

save(fullfile(golden_dir, 'unit_matlab.mat'), '-v7', '-struct', 'u');

% A small MATLAB v7.3 (HDF5) reference and trait file, to test the h5py reader
v73 = load('ref_1kgPhase3eur_LDr2p1_DEMO.mat');
keep = 1:3000;
v73.LDmat = v73.LDmat(keep, keep);
for name = {'chrnumvec', 'posvec', 'mafvec', 'is_intergenic', 'is_ambiguous'}
    v73.(name{1}) = v73.(name{1})(keep);
end
save(fullfile(golden_dir, 'ref_v73.mat'), '-v7.3', '-struct', 'v73');
save(fullfile(golden_dir, 'ref_v5.mat'), '-v7', '-struct', 'v73');
trait = load('CTG_COG_2018_DEMO.mat');
trait.logpvec = trait.logpvec(keep); trait.zvec = trait.zvec(keep);
save(fullfile(golden_dir, 'trait_v73.mat'), '-v7.3', '-struct', 'trait');
fprintf('Golden fixtures written to %s\n', golden_dir);
