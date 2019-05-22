function [logpmat,zmat]=load_gwas(traitfiles, nsnps)
%% LOAD_GWAS   Load logpvec* and zvec* from MAT files

if iscell(traitfiles), nfiles = length(traitfiles); else nfiles = 1; end;
for i=1:nfiles
    
    % Load file
    if iscell(traitfiles), file = traitfiles{i};
    else                   file = traitfiles;
    end
    fprintf('\n   %s... ',file)
    
    % Get all field names
    traits  = load(file);
    if isempty(traits),
        error( sprintf('Data file %s does not exist\n',file) );
    end
    fields = fieldnames(traits);
    
    % Find which contains single log p
    id = find(strncmpi(fields,'logp',4));
    lpfield = fields( id );
    if isempty(lpfield),
        error( sprintf('Data file %s must contain variable logp*\n',file) );
    end
    fprintf('%s... ',lpfield{1})
    lp = traits.(lpfield{1});
    
    % Find which contains z (or beta)
    id = find(strncmpi(fields,'z',1));
    zfield = fields( id );
    if isempty(zfield),
        error( sprintf('Data file %s must contain variable z*\n',file) );
    end
    fprintf('%s... ',zfield{1})
    z = traits.(zfield{1});

    if ~isvector(lp) || ~isvector(z)
        error( sprintf('logp and z must be vectors, not matrices') );
    end
    
    if length(lp) ~= nsnps, error('%s logpvec has %i snps - expected %i snps', file, length(lp), nsnps); end
    if length(z) ~= nsnps,  error('%s zvec has %i snps - expected %i snps', file, length(z), nsnps); end
    
    % Assign to logpmat 
    logpmat(:,i) = lp;
    zmat(:,i) = z; 
end


