function im_sm = SparseSmooth2d(im,wim,lamvec)
% removed colvec

[nr, nc] = size(im);

D = diff(speye(nr),2);
ivec = find(D~=0); [ivec1 ivec2] = ind2sub(size(D),ivec); vvec = D(ivec);
%kvec = colvec(repmat([1:nc],[length(ivec1) 1]));
kvec = repmat(1:nc,[length(ivec1) 1]); kvec = kvec(:);
submat = [repmat([ivec1 ivec2],[nc 1]) kvec]; valvec = repmat(vvec,[nc 1]);

L1 = spalloc((nr-2)*nc,nr*nc,nr*nc*3);
ivec1 = sub2ind([(nr-2) nc],submat(:,1),submat(:,3));
ivec2 = sub2ind([nr nc],submat(:,2),submat(:,3));
ivec = sub2ind(size(L1),ivec1,ivec2);
L1(ivec) = valvec;

D = diff(speye(nc),2);
ivec = find(D~=0); [ivec1 ivec2] = ind2sub(size(D),ivec); vvec = D(ivec);
% kvec contains repeated nc indices
% kvec = colvec(repmat([1:nr],[length(ivec1) 1]));
kvec = repmat(1:nr,[length(ivec1) 1]); kvec = kvec(:);
submat = [repmat([ivec1 ivec2],[nr 1]) kvec]; valvec = repmat(vvec,[nr 1]);

L2 = spalloc((nc-2)*nr,nr*nc,nr*nc*3);
% ivec1 and ivec2 are converted to linear indices
ivec1 = sub2ind([(nc-2) nr],submat(:,1),submat(:,3));
ivec2 = sub2ind([nr nc],submat(:,3),submat(:,2));
ivec = sub2ind(size(L2),ivec1,ivec2);
L2(ivec) = valvec;

L = cat( 1, sqrt(lamvec(1))*L1, sqrt(lamvec(1))*L2 );

wvec = wim(:);
yvec_tmp = im(:);
yvec_sm = zeros(size(yvec_tmp));

LtL = L'*L;
H = spdiags(max(0.001,wim(:)),0,nr*nc,nr*nc) + LtL;
g = wvec.*(yvec_sm-yvec_tmp) + LtL*yvec_sm;
dyvec_sm = -H\g;
yvec_sm = yvec_sm+dyvec_sm;

im_sm = reshape(yvec_sm,[nr nc]);

