function [handles, qqmat] = plot_qq_annot(logpvec1, annomat, traitname1, annonames, opts, LDmat, ...
    pruneidx, excludevec)
% 19.06 added handles as output
%       changed figure(8000 + iteri*100 + 2) to figure(8000 + iteri*100)
% 20.06 use results.pruneidx if available
% 23.06 removed colvec, rowvec
% 28.10 changed hist to histc for correctness (line 43, 55)
% 04.11 added confidence interval qqmat qqci_max qqci_min
% 09.16 added option to flip traits (i.e. plot QQ of all secondary traits conditioned on the same
% "primary" trait)

qqmat = [];

if exist('excludevec', 'var') % AMD: exclude exclude_from_fit SNPs from QQ plots
    logpvec1(excludevec) = NaN;
    annomat(excludevec, :) = NaN;
end

annothresh = 1;  % threshold for values in annomat
numannots = length(annonames);
handles = nan;
hv = opts.t1breaks;
fontsize_axes = 16;
fontsize_legends = 16;
fontsize_title = 16;

% Plotting Q-Q for each trait
scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/1.5 scrsz(4)/1.5]);
set(gcf, 'Name', sprintf('QQ %s | %s', traitname1, sprintf('%s ',annonames{:})))

h1 = nan(numannots+1,1);
for iteri = 1:numannots
    
    annovec = annomat(:, iteri);
    qqmat = zeros(length(hv), 1);
    qqci_max = zeros(length(hv), 1);
    qqci_min = zeros(length(hv), 1);
    
    % Pruning
    cntmat = zeros(length(hv),1);
    for iterj=1:size(pruneidx, 2)
        fprintf(1,'%04d/%04d... ',iterj,size(pruneidx, 2))
        
        lp1_tmp = logpvec1; lp2_tmp = annovec;

        lp1_tmp(~pruneidx(:, iterj)) = NaN;
        lp2_tmp(~pruneidx(:, iterj)) = NaN;
        for j = 1:1
            iivec = lp2_tmp >= annothresh;
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

    color1 = get(gca,'ColorOrder');
    hold all

    if iteri==1
        h1(iteri) = plot(-log10(qqmat(:,j)), hv, 'LineWidth', 2, 'color', 'k');
    else
        h1(iteri) = plot(-log10(qqmat(:,j)), hv, 'LineWidth', 2, 'color', color1(iteri-1,:));
    end
end
   
h1(numannots+1) = plot(hv, hv, 'k--');
xlim([0 7.3]); ylim([0 7.3]);
set(gca,'FontSize',fontsize_axes);
handles(1) = gcf
box(gca, 'on')

title(sprintf('SNP enrichments in %s', traitname1), 'FontSize', fontsize_title)

% Legend
legends = cell(numannots+1, 1);

for j = 1:numannots, legends{j} = annonames{j}; end
legends{numannots+1} = 'Expected';
legend(h1,legends,'Location','SouthEast','FontSize',fontsize_legends);

set(ylabel(sprintf('Nominal -log_{10}(p_{%s})',traitname1)),'FontSize',24)
set(xlabel(sprintf('Empirical -log_{10}(q_{%s})',traitname1)),'FontSize',24)

handles(iteri) = gcf;


