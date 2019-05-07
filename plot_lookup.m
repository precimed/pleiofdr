function handles = plot_lookup(logpvec1, logpmat2, traitname1, traitnames, opts, LDmat, ...
    pruneidx, excludevec, lookup12, lookup21)
% 20.06 use results.pruneidx if available
% 03.07 use lookup table if available instead of re-computing (lines 22-44)

logpthreshvec = opts.t2breaks;
hv = opts.t1breaks;
hv = hv(1:opts.thin:end);
fontsize_axes = 22;
ncondtraits = length(traitnames);
handles = nan( ncondtraits, 1 );

if opts.randprune
    if ~exist('pruneidx','var'), pruneidx=[]; end
    if isempty(pruneidx)
        % AMD: make sure we only consider SNPs with defined values when pruning
        defvec = ~excludevec & isfinite(logpvec1+sum(logpmat2,2));
        pruneidx = random_prune_idx_amd(opts.randprune_n, LDmat, defvec);
    end
end

for iteri = 1:ncondtraits
    
    logpvec2 = logpmat2(:, iteri);
    switch opts.stattype
        case 'condfdr'
            if exist('lookup12','var')
                lookup = lookup12{iteri}{1};
            else % create lookup table
                if opts.randprune
                    fprintf('\n   Trait %d/%d... ',iteri,ncondtraits)
                    lookup = lookup_table(logpvec1, logpvec2, opts, pruneidx);    
                else
                    lookup = lookup_table(logpvec1, logpvec2, opts);
                end
            end
        case 'conjfdr'
            if exist('lookup12','var') && exist('lookup21','var')
                look12 = lookup12{iteri}{1};
                look21 = lookup21{iteri}{1};
                [M,N] = size(look12);
                look12 = [ look12 ; repmat(look12(end,:), [N-M, 1]) ];
                look21 = [ look21 ; repmat(look21(end,:), [N-M, 1]) ];
                lookup = max( look12,look21' );
            else % create lookup table
                if opts.randprune
                    fprintf('\n   Trait %d/%d... ',iteri,ncondtraits)
                    lookup = conj_lookup_table(logpvec1, logpvec2, opts, pruneidx);
                else
                    lookup = conj_lookup_table(logpvec1, logpvec2, opts);
                end
            end
        otherwise
            error ('Not implemented FDR type: %s', opts.stattype);
    end
    %expstr1 = rowvec([repmat('_',[1 length(traitname1)]); traitname1]);
    %expstr2 = rowvec([repmat('_',[1 length(traitnames{iteri})]); traitnames{iteri}]);
    expstr1 = [repmat('_',[1 length(traitname1)]); traitname1];
    expstr2 = [repmat('_',[1 length(traitnames{iteri})]); traitnames{iteri}];
    expstr1 = expstr1(:)';
    expstr2 = expstr2(:)';
    
    % make figure
    figure('visible', 'on'); gcf;
    set(gcf,'Name',sprintf('Lookup %s %s',traitname1,traitnames{iteri}))
    handles(iteri) = gcf;
    if strcmpi(opts.stattype, 'condfdr')
        imagesc(lookup(:, hv<7.3), [0 1]); axis xy;
        %imagesc(-log10(lookup(:, :))); axis xy;
        colormap(flipud(hot));  
        %colormap(flipud(jet));
        colorbar
        h = title(sprintf('Conditional FDR%s _| %s',expstr1,expstr2));
        set(h,'FontSize',fontsize_axes);
        set(gca,'FontSize',fontsize_axes);
        xticks = opts.t1_low : 1 : opts.t1_up;
        set(gca,'XTick',round(interp1(hv, 1:length(hv), xticks)),'XTickLabel',xticks);
        xlabel(sprintf('-log_1_0(p%s)',expstr1),'FontSize',fontsize_axes);
        yticks = opts.t2_low : 1 : opts.t2_up;
        set(gca,'YTick',round(interp1(logpthreshvec,1:length(logpthreshvec),yticks)), ...
            'YTickLabel',yticks);
        ylabel(sprintf('-log_1_0(p%s)',expstr2),'FontSize',fontsize_axes);
    elseif strcmpi(opts.stattype, 'conjfdr')
        imagesc(lookup(hv < 7.3, hv < 7.3)', [0 1]); axis xy equal tight;
        %colormap(hot);  
        colormap(flipud(hot));
        colorbar
        h = title(sprintf('Conjunctional FDR%s _& %s',expstr1,expstr2));
        set(h,'FontSize',fontsize_axes);
        set(gca,'FontSize',fontsize_axes);
        xticks = 1:1:7;
        set(gca,'XTick',round(interp1(hv, 1:length(hv), xticks)),'XTickLabel',xticks);
        set(gca,'FontSize',fontsize_axes);
        h=xlabel(sprintf('-log_1_0(p%s)',expstr1));
        set(h,'FontSize',fontsize_axes);
        yticks = 1:1:7;
        set(gca,'YTick',round(interp1(hv, 1:length(hv), yticks)),'YTickLabel',yticks);
        set(gca,'FontSize',fontsize_axes);
        h=ylabel(sprintf('-log_1_0(p%s)',expstr2));
        set(h,'FontSize',fontsize_axes);
    else
        error ('Not implemented FDR type: %s', opts.stattype);
    end
end
