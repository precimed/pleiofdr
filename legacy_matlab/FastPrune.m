function logpvec_pruned = FastPrune(logpvec, LDmat)

% FastPrune     Prune SNPs based on 1KG LD structure.
% logpvec,      -log10(P) vector
% LDmat,        1KG LD matrix (in sparse binary format)
%
% return pruned version of logpvec with pruned SNPs set to nan.
%
% Note:
% If length of logpvec is not equal row number of LDmat you get WRONG results.
% The SNPs with largest -log10(P) will be kept in each block.
%
% logpvec_pruned = FastPrune(logpvec,LDmat);

[sv, si] = sort(abs(logpvec),'descend');
si = si(isfinite(sv));
prunevec = false(size(logpvec));

fprintf(' ')
for i = 1:length(si)

    if i==1 || mod(i, ceil(length(si)/10))==0
        out = sprintf('%d/%d.. ',floor(i),floor(length(si)));
        fprintf(out);
        for j = 1:length(out), fprintf('\b'), end
    end

    if ~prunevec(si(i))
        prunevec(LDmat(:, si(i)) ~= 0) = 1;
        prunevec(si(i)) = 0;
    end

end
logpvec_pruned = logpvec;
logpvec_pruned(prunevec) = NaN;
fprintf('\b')
