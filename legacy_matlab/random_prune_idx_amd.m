function snpidx=random_prune_idx_amd(niter, LDmat, defvec)
% 02.07 added warning for improper use of LDmat


nsnp = size(LDmat, 1);
snpidx = false(nsnp, niter);

%tic
%LDmat(1:nsnp+1:end) = 0; % to prevent "self-pruning" % AMD: what is this??

fprintf('\n   Generating random prune indices... ')
niterstrlen = length( num2str( niter ) );
fprintf( '%*d/%*d', niterstrlen, 0, niterstrlen, niter )
for k=1:niter
    fprintf('\b');
    for n=1:niterstrlen, fprintf('\b\b'); end
    fprintf( '%*d/%*d', niterstrlen, k, niterstrlen, niter )
    tmplogp = unifrnd(0, 1, [nsnp 1]);

    tmplogp(~defvec) = NaN; % AMD: remove SNPs not defined in all traits 
    
    tmplogp_p = FastPrune(tmplogp, LDmat);

    snpidx(:, k) = isfinite(tmplogp_p);
end

fprintf('\n');

%toc
end
