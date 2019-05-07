% (c) 2019 University of Oslo
% Pleiotropy-informed conditional and conjunctional false discovery rate

if ~exist('config', 'var'), config='config.txt'; end

cfg = TextConfig();

% define variable names
CFG_TRAITFOLDER_STR='traitfolder';
CFG_TRAITFILE1_STR='traitfile1';
CFG_TRAITNAME1_STR='traitname1'; 
CFG_TRAITFILES_CELL='traitfiles';
CFG_TRAITNAMES_CELL='traitnames';  
CFG_REFFILE_STR='reffile';
CFG_MATLIBRARY_PATH_STR='mlibrary';
CFG_RANDPRUNE_BOOL='randprune';
CFG_RANDPRUNE_GC_BOOL='randprune_gc';
CFG_RESET_PRUNEIDX_BOOL='reset_pruneidx';
CFG_RANDPRUNE_N_NUM='randprune_n';
CFG_RANDPRUNE_FILE_STR='randprune_file';
CFG_RANDPRUNE_REPEATS_STR='randprune_repeats';
CFG_STATTYPE_STR='stattype';
CFG_FDRTHRESH_NUM='fdrthresh';
CFG_PTHRESH_NUM='pthresh';
CFG_ONSCREEN_BOOL='onscreen';
CFG_OUTPUTDIR_STR='outputdir';
CFG_MANH_FONTSIZE_GENENAMES_NUM='manh_fontsize_genenames';
CFG_MANH_LEGEND_LOC_STR='manh_legend';
CFG_MANH_PLOT_BOOL='manh_plot';
CFG_MANH_YSPACE_NUM='manh_yspace';
CFG_MANH_YMARGIN_NUM='manh_ymargin';
CFG_MANCOLORLIST_MAT='manh_colorlist';
CFG_EXCLUDE_CHR_POS_MAT='exclude_chr_pos';
CFG_EXCLUDE_FROM_DISCOVERY_BOOL='exclude_from_discovery';
CFG_MAFTHRESH_NUM='mafthresh';
CFG_USE_STANDARD_GC_BOOL='use_standard_gc';
CFG_PERFORM_GC_BOOL='perform_gc';
CFG_EXCLUDE_AMBIGUOUS_SNPS_BOOL='exclude_ambiguous_snps';

% declare default parameters
cfg.declare(CFG_TRAITFOLDER_STR,  '../example_data_for_pleiotropy');
cfg.declare(CFG_TRAITFILE1_STR, 'PGC2_SCZ.mat');
cfg.declare(CFG_TRAITNAME1_STR, 'SCZ');
cfg.declare(CFG_TRAITFILES_CELL, '{''COG_charge.mat''}');
cfg.declare(CFG_TRAITNAMES_CELL, '{''COGNITION''}');
cfg.declare(CFG_REFFILE_STR, 'ref9545380_1kgPhase3eur_LDr2p1.mat');
cfg.declare(CFG_MATLIBRARY_PATH_STR, './');
cfg.declare(CFG_RANDPRUNE_BOOL, true);
cfg.declare(CFG_RANDPRUNE_GC_BOOL, false);
cfg.declare(CFG_RESET_PRUNEIDX_BOOL, true);
cfg.declare(CFG_RANDPRUNE_N_NUM, 20);
cfg.declare(CFG_RANDPRUNE_FILE_STR, '');
cfg.declare(CFG_RANDPRUNE_REPEATS_STR, 'default');
cfg.declare(CFG_STATTYPE_STR, 'conjfdr');
cfg.declare(CFG_FDRTHRESH_NUM, 0.05);
cfg.declare(CFG_PTHRESH_NUM, 1);
cfg.declare(CFG_ONSCREEN_BOOL, false);
cfg.declare(CFG_OUTPUTDIR_STR, 'test');
cfg.declare(CFG_MANH_FONTSIZE_GENENAMES_NUM, 12);
cfg.declare(CFG_MANH_LEGEND_LOC_STR, 'NorthEast');
cfg.declare(CFG_MANH_PLOT_BOOL, true);
cfg.declare(CFG_MANH_YSPACE_NUM, 0.75);
cfg.declare(CFG_MANH_YMARGIN_NUM, 0.25);
cfg.declare(CFG_MANCOLORLIST_MAT, ...
    '[1 0 0; 1 0.5 0 ; 0 0.75 0.75; 0 0.5 0; 0.75 0 0.75; 0 0 1; 0 1 0; 0 1 1]');
cfg.declare(CFG_EXCLUDE_CHR_POS_MAT, '[6 25119106 33854733]');  % MHC coordinates
cfg.declare(CFG_EXCLUDE_FROM_DISCOVERY_BOOL, false);
cfg.declare(CFG_MAFTHRESH_NUM, 0.005);
cfg.declare(CFG_USE_STANDARD_GC_BOOL, false);
cfg.declare(CFG_PERFORM_GC_BOOL, true);
cfg.declare(CFG_EXCLUDE_AMBIGUOUS_SNPS_BOOL, false);

% load config file if it was created
if exist(config, 'file') == 2
   cfg.loadFile(config);
end
% show current config
cfg.print();

mlibrary = cfg.get_str(CFG_MATLIBRARY_PATH_STR);  % folder containing the matlab functions
traitfolder = cfg.get_str(CFG_TRAITFOLDER_STR);
traitfile1 = cfg.get_str(CFG_TRAITFILE1_STR);
traitname1 = cfg.get_str(CFG_TRAITNAME1_STR);
traitfiles = cfg.get_cell(CFG_TRAITFILES_CELL);
traitnames = cfg.get_cell(CFG_TRAITNAMES_CELL);

addpath( mlibrary );

if ~isempty(traitfolder)
    traitfile1 = fullfile(traitfolder, traitfile1);
    for i=1:length(traitfiles), traitfiles{i} = fullfile(traitfolder, traitfiles{i}); end
end

clear options
options = pleioOpt();
options.randprune = cfg.get_bool(CFG_RANDPRUNE_BOOL);
options.randprune_gc = cfg.get_bool(CFG_RANDPRUNE_GC_BOOL);
options.reset_pruneidx = cfg.get_bool(CFG_RESET_PRUNEIDX_BOOL);
options.randprune_n = cfg.get_num(CFG_RANDPRUNE_N_NUM);
options.randprune_file = cfg.get_str(CFG_RANDPRUNE_FILE_STR);
options.randprune_repeats = cfg.get_str(CFG_RANDPRUNE_REPEATS_STR);
options.stattype = cfg.get_str(CFG_STATTYPE_STR);
options.fdrthresh = cfg.get_num(CFG_FDRTHRESH_NUM);
options.pthresh = cfg.get_num(CFG_PTHRESH_NUM);
options.onscreen = cfg.get_bool(CFG_ONSCREEN_BOOL);
options.outputdir = cfg.get_str(CFG_OUTPUTDIR_STR);
options.manh_fontsize_genenames = cfg.get_num(CFG_MANH_FONTSIZE_GENENAMES_NUM);
options.manh_plot = cfg.get_bool(CFG_MANH_PLOT_BOOL);
options.manh_legend = cfg.get_str(CFG_MANH_LEGEND_LOC_STR);
options.manh_yspace = cfg.get_num(CFG_MANH_YSPACE_NUM);
options.manh_ymargin = cfg.get_num(CFG_MANH_YMARGIN_NUM);
options.manh_colorlist = 0.8 * cfg.get_mat(CFG_MANCOLORLIST_MAT);
options.exclude_from_fit_and_discovery = cfg.get_bool(CFG_EXCLUDE_FROM_DISCOVERY_BOOL);
options.use_standard_gc = cfg.get_bool(CFG_USE_STANDARD_GC_BOOL);
options.perform_gc = cfg.get_bool(CFG_PERFORM_GC_BOOL);
options.exclude_ambiguous_snps = cfg.get_bool(CFG_EXCLUDE_AMBIGUOUS_SNPS_BOOL);
options.mafthresh = cfg.get_num(CFG_MAFTHRESH_NUM);

%% LOAD FILES

if ~exist('LDmat', 'var')
    reffile=cfg.get_str(CFG_REFFILE_STR);
    fprintf('Loading reference file ("%s") ...', reffile);
    load(reffile);
    % make sure everything went well
    assert(exist('LDmat', 'var')  == 1,  'error loading LD-matrix');
    assert(ismatrix(LDmat),              'error loading LD-matrix');
    assert(exist('mafvec', 'var') == 1,  'error loading mafvec');
    assert(ismatrix(mafvec),             'error loading mafvec');
    assert(exist('is_intergenic', 'var') == 1, 'error loading is_intergenic');
    ivec0=logical(is_intergenic);
    fprintf('done\n')
end

if ~strcmp( options.randprune_file, '' ) && exist( options.randprune_file, 'file' )
    fprintf('Loading prune set %s... ', options.randprune_file);
    loaded_data = load( options.randprune_file );
    tmpfields = fieldnames( loaded_data );
    pruneidx = loaded_data.( tmpfields{1} );
    options.randprune_n = size( pruneidx, 2 );
    options.reset_pruneidx = false;
    fprintf('done\n')
end

exclude_chr_pos=cfg.get_mat(CFG_EXCLUDE_CHR_POS_MAT);
for i=1:size(exclude_chr_pos, 1)
    from_idx = find(chrnumvec == exclude_chr_pos(i, 1) & posvec >= exclude_chr_pos(i, 2), 1 );
    to_idx = find(chrnumvec == exclude_chr_pos(i, 1) & posvec <= exclude_chr_pos(i, 3), 1, 'last' );
    options.exclude_from_fit{end+1} = [from_idx, to_idx];
end

if ~options.onscreen, set(0, 'DefaultFigureVisible', 'off'); end

%% Execute pleiotropy analysis, make all figures and resulting tables
pleiotropy_analysis
