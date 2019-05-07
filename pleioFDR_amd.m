function [fdrmat, fdrvec0, pruneidx, fdrmat12, fdrmat21, lookup12, lookup21] =  ...
   pleioFDR_amd(logpvec1, logpmat2, opts, LDmat, excludevec, pruneidx)
%% pleioFDR  Perform pleiotropy informed False Discovery Rate method
%     using conditional or conjuction FDR. Calls on cond_FDR to calculate 
%     phenotype 1 conditioned on phenotype(s) 2. Conjuction FDR 
%     additionally calls cond_FDR for phenotype(s) 2 on 1 and
%     finds the max. 
%  Inputs:
%     logpvec1   : (#SNP x 1) -log10(p-value) of phenotype 1
%     logpmat2   : (#SNP x M) -log10(p-value) of M phenotype(s)
%                  to condition on
%     opts       : Options (see pleioOpt.m)
%                  Select opts.stattype as 'condfdr' or 'conjfdr'
%     LDmat      : (#SNP x #SNP) linkage-disequilibrium matrix
%     excludevec : Optional (#SNP x 1) binary list of SNPs to exclude
%  Outputs:
%     fdrmat     : (#SNP x M) FDR matrix conditioned/conjuctioned on M phenotypes 
%     fdrvec0    : (#SNP x 1) unconditioned FDR of phenotype 1
%     fdrmat12, fdrmat21  : extra outputs for conjunction FDR (1 on 2, 2 on 1)
%  Required functions :
%     random_prune_idx.m, cond_FDR.m
%  References :
%     Andreassen, O. A., Djurovic, S., Thompson, W. K., Schork,
%     A. J., Kendler, K. S., O?Donovan, M. C., ... & Dale, A. M. (2013).
%     Improved detection of common variants associated with schizophrenia
%     by leveraging pleiotropy with cardiovascular-disease risk factors.
%     The American Journal of Human Genetics, 92(2), 197-209.
%     http://www.ncbi.nlm.nih.gov/pubmed/23375658
%
%  02.07.14: changes to output lookup
%  06.08.14: removes checks to exclude_from_fit for easier use
%  04.11.14: take pruneidx as input
%  11.11.15: creates pruneidx if input pruneidx is empty (line 46) 
%  18.11.15: fixed bug of pruneidx

% AMD: Move check for excludevec into cond_FDR_amd, to exclude from creation of FDR lookup table,
% not FDR value assignment
if exist('excludevec', 'var')
%    logpvec1(excludevec) = NaN;
%    logpmat2(excludevec, :) = NaN;
end
if ~exist('pruneidx','var')
    pruneidx = [];
end

ncondtraits = size(logpmat2, 2);
fdrmat12 = NaN(size(logpvec1,1), ncondtraits); fdrmat21 = fdrmat12;
lookup12 = []; lookup21 = [];
if opts.randprune
    
    if isempty(pruneidx)
        % AMD: make sure we only consider SNPs with defined values when pruning
        defvec = ~excludevec & isfinite(logpvec1+sum(logpmat2,2));
        pruneidx = random_prune_idx_amd(opts.randprune_n, LDmat, defvec);
    end
    
    for iteri=1:ncondtraits
        fprintf('\n   Trait %d/%d... ',iteri,ncondtraits)
        [ fdrmat12(:, iteri), lookup12{iteri}, fdrmat0(:,iteri) ] = ...
            cond_FDR_amd(logpvec1, logpmat2(:,iteri), opts, pruneidx, excludevec);
        if strcmpi(opts.stattype, 'conjfdr')
            [ fdrmat21(:, iteri), lookup21{iteri}, fdrmat02(:,iteri) ] = ...
                cond_FDR_amd(logpmat2(:,iteri), logpvec1, opts, pruneidx, excludevec);
        end
    end
    %fdrvec0_pruned = mean(fdrmat0,2); 
    %if strcmpi(opts.stattype, 'conjfdr')
    %    fdrvec02_pruned = mean(fdrmat02,2); 
    %end
    
    % the unconditioned fdr is calculated from unpruned SNPs    
    [~ , ~, fdrvec0] = cond_FDR_amd(logpvec1, logpmat2(:,1), opts, [], []);
    
else
    
    pruneidx = [];
    for iteri=1:ncondtraits
        [fdrmat12(:,iteri), lookup12{iteri}, fdrvec0] = ...
            cond_FDR_amd(logpvec1, logpmat2(:, iteri), opts, pruneidx, excludevec);
        if strcmpi(opts.stattype, 'conjfdr')
            [fdrmat21(:,iteri), lookup21{iteri}, fdrvec02] = ...
                cond_FDR_amd(logpmat2(:,iteri), logpvec1, opts, pruneidx, excludevec);
        end
    end
    
end

switch opts.stattype
    case 'condfdr'
        fdrmat = fdrmat12;
    case 'conjfdr'
        fdrmat = max(fdrmat12, fdrmat21);
    otherwise
    error ('Not-implemented FDR type %s', opts.stattype);
end
