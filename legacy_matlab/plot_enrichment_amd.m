function [handles, enrichmat] = plot_enrichment_amd(logpvec1, logpmat2, traitname1, traitnames, opts, LDmat, ...
    pruneidx, excludevec, flip_traits)
% 19.06 added handles as output
% 23.06 removed colvec, rowvec
% AMD: added excludevec explicitly in Q-Q plots
% 09.16 added option to flip traits (i.e. plot QQ of all secondary traits conditioned on the same
% "primary" trait)

enrichmat = [];

ncondtraits = length(traitnames);
handles = nan( ncondtraits, 1 );
logpthreshvec = opts.qqbreaks;
hv = opts.t1breaks;
fontsize_legends = 16;
fontsize_title = 16;

if opts.randprune
    if ~exist('pruneidx','var'), pruneidx=[]; end
    if isempty(pruneidx)
        % AMD: make sure we only consider SNPs with defined values when pruning
        defvec = ~excludevec & isfinite(logpvec1+sum(logpmat2,2));
        pruneidx = random_prune_idx_amd(opts.randprune_n, LDmat, defvec);
    end
end

if exist('excludevec', 'var') % AMD: exclude exclude_from_fit SNPs from QQ plots
    logpvec1(excludevec) = NaN;
    logpmat2(excludevec, :) = NaN;
end
if ~exist('flip_traits', 'var')
    flip_traits = false;
end

% Plotting enrichment for each trait
scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/1.5 scrsz(4)/1.5]);
if ~flip_traits
    set(gcf, 'Name', sprintf('Enrichment %s | %s', traitname1, sprintf('%s ', traitnames{:})))
else
    set(gcf, 'Name', sprintf('Enrichment %s | %s', sprintf('%s ', traitnames{:}), traitname1))
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
                
                %chc = cumsum(hc)/sum(hc);
                %tmp = 1 - chc(:);
                %defvec = isfinite(tmp);
                %qqmat(defvec, j) = qqmat(defvec, j) + tmp(defvec);
                %cntmat(:,j) = cntmat(:, j) + double(defvec);
                
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
            [ phat, pci ] = binofit_wrap(cumsum(hc), sum(hc));
            qqmat(:,j) = 1-phat(:);
            qqci_max(:,j) = 1-pci(:,1);
            qqci_min(:,j) = 1-pci(:,2);
            
        end
    end
    if ~flip_traits
        trait1 = traitname1; trait2 = traitnames{iteri};
    else
        trait1 = traitnames{iteri}; trait2 = traitname1;
    end

    % Legend
    legends = cell( length(logpthreshvec), 1 );
    for j = 1:length(logpthreshvec)
        if logpthreshvec(j) == 0
            legends{j} = 'All SNPs';
        else
            legends{j} = sprintf('p_{%s} < 10^{-%d}', trait2, logpthreshvec(j));
            
        end
    end

%   enrichmat = bsxfun(@minus, -log10(qqmat(:,1)), -log10(qqmat));
    enrichmat = bsxfun(@rdivide, qqmat, qqmat(:,1) );

    if opts.show_ci
%         enrichmax = bsxfun(@minus, -log10(qqmat(:,1)), -log10(qqci_max));
%         enrichmin = bsxfun(@minus, -log10(qqmat(:,1)), -log10(qqci_min));
        enrichmax = bsxfun(@rdivide, qqci_max, qqmat(:,1) );
        enrichmin = bsxfun(@rdivide, qqci_min, qqmat(:,1) );
    end
    
    % Plot enrichment as subplots
    spcols = ceil(sqrt(ncondtraits));
    sprows = ceil(ncondtraits/spcols);
    subplot(sprows,spcols,iteri)
    
    %plot(hv, 10 .^ enrichmat,'LineWidth',2);
    color1 = get(gca,'ColorOrder');
    hold off
    htmp = nan(size(qqmat,2), 1);
    for j=1:size(qqmat,2)
        htmp(j) = plot(hv,enrichmat(:,j),'LineWidth',2,'color',color1(j,:));
        hold all
        if opts.show_ci
            plot(hv,enrichmax(:,j),'--','LineWidth',1,'color',color1(j,:));
            plot(hv,enrichmin(:,j),'--','LineWidth',1,'color',color1(j,:));
        end
    end
    xlim([0 7.3]);
    if (~is_octave())
        h=legend(htmp,legends,'Location','NorthWest');
        set(h,'FontSize',fontsize_legends);
    end
    title(sprintf('%s | %s',trait1,trait2),'FontSize',fontsize_title)
    handles(iteri) = gcf;
    if (flip_traits)
        set(ylabel(sprintf('Fold Enrichment %s | %s',trait1,trait2)),'FontSize',24)
        set(xlabel(sprintf('Nominal -log_{10}(p_{%s})',trait1)),'FontSize',24)
    else
        set(ylabel(sprintf('Fold Enrichment %s | %s',trait2,trait1)),'FontSize',24)
        set(xlabel(sprintf('Nominal -log_{10}(p_{%s})',trait2)),'FontSize',24)
    end
    
end

handles = handles( isfinite( handles ) );

% adapting maximum y
ymax = 0;
for iteri=1:ncondtraits
    subplot(sprows,spcols,iteri)
    ymax = max(ymax, max(ylim));
end
for iteri=1:ncondtraits
    subplot(sprows,spcols,iteri)
    ylim([0,ymax]);
    set(gca,'YScale','linear')
%    set(gca,'YScale','log')
end
