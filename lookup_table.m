function [looktable, lookcount, lookstd] = lookup_table(lp1, lp2, opts, pruneidx)
% the conditional FDR is estimated via conditioning on the second phenotype
% via a grid / lookup table
% 21.06 : recode lookup table using hist3 for efficiency
%         logit/sigmoid as inline function, remove dependency of logit.m
%         rewrite tmp_lp1 and tmp_lp2 (line 24) for efficiency
% 02.07 : added outputs lookcount and lookstd
% 19.08 : t1breaks and t2breaks defined as cell in pleioOpt.m
%         changed to t1breaks{1} and t2breaks{1} at line 74 and line 115
% 19.08:  changed the muim trim, trim at x and n instead (we can trim the
%         cumulative histogram)
% 10.10 : added options.smooth_lookup to turn on/off smoothing
%         (change to pleioOpt)
% 05.01 : edited options.adjust_lookup for adjusted p-values (to ensure
%          monotonicity of FDR lookup tables)




if exist('pruneidx', 'var')

    if size(pruneidx, 1) ~= size(lp1, 1)
        error('Random prune index does not conform to logpvec 1');
    end
    if size(pruneidx, 1) ~= size(lp2, 1)
        error('Random prune index does not conform to logpvec 2');
    end
    
    
    looktable = zeros(opts.t2_nbreaks, length(opts.t1breaks(1:opts.thin:end)));
    lookcount = looktable;
    lookvar   = looktable;
    fprintf('Fill lookup table 000/%3d ', opts.randprune_n);
    for iteri=1:opts.randprune_n
        fprintf('\b\b\b\b\b\b\b\b')
        fprintf('%3d/%3d ',iteri, opts.randprune_n);
        prunei = pruneidx(:,iteri);
        tmp_lp1 = NaN(size(lp1)); tmp_lp1(prunei) = lp1(prunei);
        tmp_lp2 = NaN(size(lp2)); tmp_lp2(prunei) = lp2(prunei);
        
        tmp_look = ind_look(tmp_lp1, tmp_lp2, opts);
        defmat = isfinite(tmp_look);
        lookcount = lookcount + double(defmat);
        looktable(defmat) = looktable(defmat) + tmp_look(defmat) ;
        lookvar(defmat) = lookvar(defmat) + tmp_look(defmat).^2 ;
    end
    %fprintf('\n')
    looktable = looktable ./ lookcount;
    lookstd = sqrt((lookvar./lookcount) - looktable.^2); % biased estimator
        
else
    
    looktable = ind_look(lp1, lp2, opts);
    lookcount = [];
    lookstd = [];
    
end


% FDR Adjusted p-values? (to ensure monotonicity)
% This needs some rethinking...
if opts.adjust_lookup
    fprintf('adjusting lookup... ')
    for cols = 1: size(looktable,2)  %:-1:1
        for rows = 2: size(looktable,1)  %:-1:1
            looktable(rows,cols) = ...
                min(looktable(rows,cols), looktable(rows-1,cols));
        end
        if cols>1
            looktable(:,cols) = min(looktable(:,cols),looktable(:,cols-1));
        end
    end
end

end



function tmp_look = ind_look(lp1, lp2, opts)

% hcmat12 = NaN(opts.t2_nbreaks, opts.t1_nbreaks);
% for j = 1:(opts.t2_nbreaks)
%
%     %iivec indexes all SNPs that have non-missing values in both traits
%     %and a p-value in the second trait ABOVE the logpthresvec(j)
%     iivec = (isfinite(lp1 + lp2) & lp2 >= opts.t2breaks(j));
%
%     %hc histogramm counts based on grid hv (length 1001)
%     %hc = hist(lp1(iivec), opts.t1breaks);
%     hc = histc(lp1(iivec), opts.t1breaks);
%     chc = cumsum(hc)/sum(hc);
%
%     %FDR
%     hcmat12(j,:) = hc;
% end

% more efficient codes to create lookup histogram table hcmat12(i,j)
% where opts.t2breaks(i) <= lp2 AND opts.t1breaks(j) <= lp1 < opts.t2breaks(j+1)
lpmat = [lp2,lp1]; lpmat( any(isnan(lpmat),2),: )=[];
edges = { [opts.t2breaks,inf], [opts.t1breaks,inf] };  % infinity to capture lp2>max(opts.t2breaks)

edges{2} = edges{2} - 0.5*(edges{2}(2)-edges{2}(1));  % shift half-bin to imitate old codes

hlp  = hist3(lpmat,'Edges',edges);
if is_octave()
    hlp(1,1002)=0;
    hlp(32,1)=0;
end
hcmat12 = flipud( cumsum(flipud(hlp),1) );  % lp2-cumulative histogram

% binomialfit for the histogram counts
% estimate p = cumulative counts of event / all counts
% cumulative (over trait1) counts of events  =  1+colvec(cumsum(hcmat12,2))
% all counts (marginal sum) = 1+colvec(repmat(sum(hcmat12,2)
% size(cumsum(hcmat12,2)) = 31 1001
% pci confidence intervalls for p

%x = cumsum(hcmat12, 2);   %<-- is this correct?
x = fliplr( cumsum(fliplr(hcmat12),2) );
n = repmat( sum(hcmat12,2), [1, size(x,2)] );

%trim to conform to legacy codes
x=x(1:end-1,1:end-1);
n=n(1:end-1,1:end-1);

[pest, pci] = binofit_wrap(1+x(:),1+n(:));
pest = reshape(pest, size(x));
pupp = reshape(pci(:,2), size(x));
plow = reshape(pci(:,1), size(x));

%smooth pest in the logit domain

%muim reshapes pest (vector) into muim (matrix), and does the logit transform
logit = @(x) log(x./(1-x));
muim = logit(pest);
if (opts.smooth_lookup)
    wim = (logit(pupp)-logit(plow)).^-2;
    %muim = muim(1:end-1,1:end-1); wim = wim(1:end-1,1:end-1); %trim to conform
    %to legacy codes : note move trim up to line 94
    muim = muim(:,1:opts.thin:end); wim = wim(:, 1:opts.thin:end);
    wim(~isfinite(muim)) = 0; muim(~isfinite(muim)) = 0;
    muim_sm = SparseSmooth2d(muim, wim, [opts.smf opts.smf]);
else
    muim_sm = muim(:,1:opts.thin:end);
    disp('WARNING!!! no smoothing')
end

% qim is the actually observed F(p_1 \mid p_2)
%qim = 1-logit(muim_sm,true);
sigmoid = @(x) 1./(1+exp(-x));
%qim = 1- sigmoid(muim_sm);   % <- see definition x line 77
qim = sigmoid(muim_sm);

% pim is the F_0 under the Null
hv = opts.t1breaks(1:opts.thin:end);
pim = repmat(10.^-hv,[size(qim,1) 1]);

tmp_look = min(1, max(0, pim ./ qim));

end
