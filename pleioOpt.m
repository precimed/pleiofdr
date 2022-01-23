classdef pleioOpt
    % Class wrapping options used by pleiotropy cond/conj analysis
    % -----------------
    % Properties:
    % ----------
    % stattype: 'condfdr'
    % fdrthresh: 0.0500
    % pthresh: 5.0000e-08
    % exclude_from_fit: [1 20]
    % onscreen: false
    % t1_low: 0
    % t1_up: 10
    % t1_nbreaks: 1001
    % t2_low: 0
    % t2_up: 3
    % t2_nbreaks: 31
    % qq_low: 2
    % qq_up: 3
    % qq_nbreaks: 4
    % qqbreaks: [2 2.3333 2.6667 3]
    % t1breaks:
    % t2breaks:
    % 
    % Note:
    % ----
    % No methods implemented for manipulating the objects
    %
    % Notes: 
    %   23.06  added options for manhattan plots
    %   30.06  added options for colorlist
    %   19.08  IMPORTANT! change t1breaks and t2breaks into cell so they can be modified after
    %          definition
    %   19.08  renamed t1breaks as trait2breaks to avoid confusion later
    %   10.10  added options.smooth_lookup to turn on/off smoothing (see lookup_table.m)
    %   15.01  added show_ci for confidence intervals in q-q and enrichment plots
    
    properties % default paramters
        stattype = 'condfdr'; 
        fdrthresh = 0.05; pthresh = 5e-8; exclude_from_fit = {}; 
        mafthresh = 0.005;
        t1_low = 0; t1_up = 30; t1_nbreaks = 3001;
        t2_low = 0; t2_up = 3;  t2_nbreaks = 31;
        qq_low = 0; qq_up = 3; qq_nbreaks = 4;
        onscreen = false;
        randprune = true;
        randprune_n = 100;
        randprune_file = '';
        randprune_repeats = 'default';  % choose from 'default', 'maxout', 'none'
        reset_pruneidx = true;
        smf = 1e2;
        thin=10;
        perform_gc = true;
        randprune_gc = false;
        save_figs_as_group = true;
        smooth_lookup = true;
        adjust_lookup = false;
        show_ci = false;         %
        qq_annotation_index = [];
        qq_annotation_r2 = [];
        use_standard_gc = false;
        exclude_ambiguous_snps = false;
        dummy_zscore = false;

        % Controls whether SNPs from exclude_from_fit regions should be
        % excluded from FDR fit only, or also from discovery
        % Set this option to TRUE will exclude SNPs from FDR table,
        % manhattan plot, etc.
        exclude_from_fit_and_discovery = false;

        outputdir = '';                         % Folder to place the result
        correct_for_sample_overlap = false;     % Boolean flag indicating whether to control for sample overlap
        fishercomb  = true;                     % Boolean flag indicating whether to use combined fisher statistics

        % manhattan plot options
        manh_fontsize_genenames = 20; %16;      % (will be rescaled later)
        manh_fontsize_legends = 18; %22 
        manh_fontsize_axes = 18; %22
        manh_legend = 'NorthEast';
        manh_plot = false;
        manh_redraw = false;                    % frequent redraw to see progress (slower)
        manh_ymargin = 1; %0.25;                % margin between points and genenames
        manh_yspace = 0.75 ; %0.5               % line spacing between genenames
        manh_image_size = [1920, 1440];         % info about image resolution used to filter
                                                % points on manhatten plot;
                                                % NaN means not to filter.

        manh_colorlist = 0.8 *...               % set lower for darker colors
         [  1.00 0    0   ;
            1.00 0.5  0   ;
            0    0.75 0.75;
            0    0.50 0   ;
            0.75 0    0.75;
            0    0    1.00 ]; 
        %    0.25 0.25 0.25;
        %    0.5  0.5  0.5 ;  
        %    0    1.0  0   ;
        %    0    0    0.5 ]
            
    end
    properties (Dependent=true)
        qqbreaks; t1breaks; t2breaks;
    end
    methods
        function obj=set.stattype(obj, stype)
            if ~(strcmpi(stype, 'condfdr') ||...
                    strcmpi(stype, 'conjfdr'))
                error('Not-implemented FDR type %s', stype);
            end
            obj.stattype = stype;
        end
        function obj=set.exclude_from_fit(obj, regionidx)
            for j=1:length(regionidx)
                if regionidx{j}(1)  >= regionidx{j}(2)
                    error('exclude region is empty');
                end
                if any(regionidx{j} <= 0)
                    error('Invalid exclude region: negative values');
                end
            end
            obj.exclude_from_fit = regionidx;
        end
        function tmp=get.qqbreaks(obj)
            tmp=linspace(obj.qq_low, obj.qq_up, obj.qq_nbreaks);
        end
        function obj=get.t1breaks(obj)
            obj=linspace(obj.t1_low, obj.t1_up, obj.t1_nbreaks);
        end
        function obj=get.t2breaks(obj)
            obj=linspace(obj.t2_low, obj.t2_up, obj.t2_nbreaks);
        end


    end
end
