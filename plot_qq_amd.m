function [handles, qqmat] = plot_qq_amd(logpvec1, logpmat2, traitname1, traitnames, opts, LDmat, ...
    pruneidx, excludevec, flip_traits, use_annotations, plot_tdr)
% 19.06 added handles as output
%       changed figure(8000 + iteri*100 + 2) to figure(8000 + iteri*100)
% 20.06 use results.pruneidx if available
% 23.06 removed colvec, rowvec
% 28.10 changed hist to histc for correctness (line 43, 55)
% 04.11 added confidence interval qqmat qqci_max qqci_min
% 09.16 added option to flip traits (i.e. plot QQ of all secondary traits conditioned on the same
% "primary" trait)

qqmat = [];

if nargin<3
    traitname1 = 'trait 1';
    traitnames = 'trait 2';
end 
if nargin<5
    opts = struct('randprune',false,'t1breaks',linspace(0,10,1001));
end

if exist('excludevec', 'var') % AMD: exclude exclude_from_fit SNPs from QQ plots
    logpvec1(excludevec) = NaN;
    logpmat2(excludevec, :) = NaN;
end

if ~exist('flip_traits', 'var'), flip_traits = false; end
if ~exist('use_annotations', 'var'), use_annotations = false; end
if ~exist('plot_tdr', 'var'), plot_tdr = false; end

ncondtraits = length(traitnames);
handles = nan( ncondtraits, 1 );
if ~use_annotations, logpthreshvec = opts.qqbreaks;
else, logpthreshvec = opts.qq_annotation_r2; end
if isempty(logpthreshvec), return; end
hv = opts.t1breaks;
fontsize_axes = 20;
fontsize_legends = 20;
fontsize_title = 26;

% Create random pruning index if not available
if opts.randprune
    if ~exist('pruneidx','var'), pruneidx=[]; end
    if isempty(pruneidx)
        % AMD: make sure we only consider SNPs with defined values when pruning
        defvec = ~excludevec & isfinite(logpvec1+sum(logpmat2,2));
        pruneidx = random_prune_idx_amd(opts.randprune_n, LDmat, defvec);
    end
end

% Plotting Q-Q for each trait
scrsz = [0 0 1920 1080]
figure('Position',[1 scrsz(4)/2 scrsz(4)/1.5 scrsz(4)/1.5]);
if ~flip_traits
    set(gcf, 'Name', sprintf('QQ %s | %s', traitname1, sprintf('%s ',traitnames{:})))
else
    set(gcf, 'Name', sprintf('QQ %s | %s', sprintf('%s ', traitnames{:}), traitname1))
end

for iteri = 1:ncondtraits
    
    logpvec2 = logpmat2(:, iteri);
    qqmat = zeros(length(hv), length(logpthreshvec));
    qqci_max = zeros(length(hv), length(logpthreshvec));
    qqci_min = zeros(length(hv), length(logpthreshvec));
    
    % Pruning
    if opts.randprune
        cntmat = zeros(length(hv),length(logpthreshvec));
        for iterj=1:opts.randprune_n
            fprintf(1,'%04d/%04d... ',iterj,opts.randprune_n)
            
            if ~flip_traits
                lp1_tmp = logpvec1; lp2_tmp = logpvec2;
            else
                lp1_tmp = logpvec2; lp2_tmp = logpvec1;
            end

            lp1_tmp(~pruneidx(:, iterj)) = NaN;
            lp2_tmp(~pruneidx(:, iterj)) = NaN;
            for j = 1:length(logpthreshvec)
                iivec = lp2_tmp >= logpthreshvec(j);
                hc = histc(lp1_tmp(iivec), hv); hc = hc(:);
                
                [ phat, pci ] = binofit_wrap(cumsum(hc), sum(hc));
                defvec = isfinite(phat(:));
                if ~any(defvec), continue; end
                qqmat(defvec, j)    = qqmat(defvec, j) + 1-phat(:);
                qqci_max(defvec, j) = qqmat(defvec, j) + 1-pci(:,1);
                qqci_min(defvec, j) = qqmat(defvec, j) + 1-pci(:,2);
                cntmat(defvec, j)   = cntmat(defvec, j) + double(defvec);
            end
            for backs=1:13, fprintf(1,'\b'), end
        end
        qqmat = qqmat ./ cntmat;
        qqci_max = qqci_max ./ cntmat;
        qqci_min = qqci_min ./ cntmat;
    else
        for j = 1:length(logpthreshvec)
            if ~flip_traits
                iivec = logpvec2 >= logpthreshvec(j);
                hc = histc(logpvec1(iivec), hv); hc = hc(:);
            else
                iivec = logpvec1 >= logpthreshvec(j);
                hc = histc(logpvec2(iivec), hv); hc = hc(:);
            end
            %qqmat(:,j) = 1 - cumsum(hc)/sum(hc);
            
            [ phat, pci ] = binofit_wrap(cumsum(hc), sum(hc));
            qqmat(:,j) = 1-phat(:);
            qqci_max(:,j) = 1-pci(:,1);
            qqci_min(:,j) = 1-pci(:,2);
        end
    end

    % Plot QQ as subplots
    spcols = ceil(sqrt(ncondtraits));
    sprows = ceil(ncondtraits/spcols);
    subplot( sprows, spcols, iteri );
    %plot(-log10(qqmat), hv, 'LineWidth', 2);
    color1 = get(gca,'ColorOrder');
    hold off
    h1 = nan(size(qqmat,2),1);

    for j=1:size(qqmat,2)
        if plot_tdr
            tdr_vec = 1-10.^(-hv) ./ [1, qqmat(1:end-1, j)'];
            tdr_vec(tdr_vec < 0) = 0;
            h1(j) = plot(hv, tdr_vec, 'LineWidth', 2, 'color', color1(j,:));
            hold all
            continue
        end
        h1(j) = plot(-log10(qqmat(:,j)), hv, 'LineWidth', 2, 'color', color1(j,:));
        hold all
        if opts.show_ci
            plot(-log10(qqci_max(:,j)), hv, '--', 'LineWidth', 1, 'color', color1(j,:));
            plot(-log10(qqci_min(:,j)), hv, '--', 'LineWidth', 1, 'color', color1(j,:));
        end
    end
    if ~plot_tdr
        h1(j+1) = plot(hv, hv, 'k--', 'LineWidth', 1.5);
        %axis equal 
        xlim([0 7.3]); ylim([0 7.3]);
    else
        xlim([0 7.3]); ylim([0 1]);
    end
    %xl = sprintf('Empirical -log_{10}(q_{%s | %s})',traitname1,traitnames{iteri});
    %xlabel(xl,'FontSize',fontsize_axes);
    %yl = sprintf('Nominal -log_{10}(p_{%s})',traitname1);
    %ylabel(yl,'FontSize',fontsize_axes);
    set(gca,'FontSize',fontsize_axes);
    if ~flip_traits
        trait1 = traitname1; trait2 = traitnames{iteri};
    else
        trait1 = traitnames{iteri}; trait2 = traitname1;
    end

    title(sprintf('%s | %s',trait1,trait2), 'fontweight', 'normal', 'FontSize', fontsize_title)

    % Legend
    legends = cell(length(logpthreshvec), 1);
    if 1 %(iteri == ncondtraits)
        for j = 1:length(logpthreshvec)
            if logpthreshvec(j) == 0
                legends{j} = 'All SNPs';
            else
                if use_annotations
                    legends{j} = sprintf('LD_{%s} > %d', trait2, logpthreshvec(j));
                else
                    legends{j} = sprintf('p_{%s} < 10^{-%d}', trait2, logpthreshvec(j));
                end
            end
        end
        legends{j+1} = 'Expected';
        if (~is_octave())
            legend(h1,legends,'Box','off','Location','SouthEast','FontSize',fontsize_legends);
        end
    end
    handles(iteri) = gcf;
    set(gca,'linewidth',1)
    set(gca,'TickLength',[0.014 0.014])
    set(gca,'XTick',[0:1:7])
    if plot_tdr
        set(xlabel(sprintf('Nominal -log_{10}(p_{%s})',trait1)),'FontSize',24)
        set(ylabel(sprintf('Conditional TDR_{%s|%s}',trait1, trait2)),'FontSize',24)
    else
        set(gca,'YTick',[0:1:7])
        set(ylabel(sprintf('Nominal -log_{10}(p_{%s})',trait1)),'FontSize',24)
        set(xlabel(sprintf('Empirical -log_{10}(q_{%s})',trait1)),'FontSize',24)
    end
end

handles = handles( isfinite( handles ) );


