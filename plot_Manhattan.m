function  handles = plot_Manhattan(results,traitname1,traitnames,chrnumvec,options)
%% PLOT_MANHATTAN Make Manhattan plots based on
%     results.fdrmat    : (SNP x trait) FDR matrix
%     results.imat      : (SNP x 1) SNPs above -log10(fdrthresh)
%     results.imat2     : (SNP x 1) SNP loci after pruning
%     options.stattype  : condfdr/conjfdr
%     options.fdrthresh : FDR threshold
%
% Remember to update that imat, imat2 before plotting
% [imat, imat2, logfdrmat] = ind_loci_idx(fdrmat, flp, LDmat, options);
%
% Developer notes
% 13.06: changed for-loop plot phase
%        removed traitsequences option, just use one traitsequence
% 16.06: major cleanup
%        removed genenamesubst
%        gnlist setorder to 'stable', renamed extmat to textbox
%        renamed logpmat_sort into logfdrmat as it should be
%        removed p-values (logpvec1, logpmat2), use only results.fdrmat 
% 19.06: added handles
% 23.06: added manhattan plot options to options
% 25.06: rescales fontsize.genenames when rescaling y-axis (line 278)
% 30.06: added colorlist as options.manh_colorlist (update pleioOpt.m)
% 02.07: fixed bug for no genes passing threshold (line 274)
% 22.05: don't plot fdrvec0 if snp is NaN on trait2
%        fdrvec0( any(isnan(fdrmat),2) ) = NaN; (line 169)


%% Customizable plot options

chrnumlist = { double(unique(chrnumvec)') };            % chromosome subsets
exclude_locations = { };            % start/end points of excluded locations, in cell of arrays
showgenes = false;                     % switch off for faster plots
fontsize_genenames = 16; %14; %16;
fontsize_legends   = 18; %18; %22
fontsize_axes      = 13; % 18; %22
fontsize_label      = 18; % 18; %22
ymargin = options.manh_ymargin; % 1; %0.25;      % margin between points and genenames
yspace  = options.manh_yspace; % 0.75; %0.5     % line spacing between genenames
LegendLocation = 'NorthEast';
graylevel = 0.1;                   % grayness level for chromosome fill
colorlist = options.manh_colorlist;

%% Legacy codes

fdrmat       = results.fdrmat;
imat         = results.imat;
imat2        = results.imat2;
chrnumvec    = results.chrnumvec;
genenamelist = results.genenamelist;
fdrvec0      = results.fdrvec0;

ntraits  = size(fdrmat, 2);
%logpmat1 = repmat(logpvec1,[1,ntraits]);
%defmat   = isfinite(logpmat1+logpmat2);
logfdrmat  = -log10(fdrmat);
%logfdrmat(~defmat) = NaN;


excludevec = false(size(logfdrmat,1),1);
for j = 1:length(exclude_locations)
    excludevec(exclude_locations{j}(1):exclude_locations{j}(2)) = true;
end


%% PREPARE LEGEND AND TRAIT PLOTTING SEQUENCE -----------------------

switch(options.stattype)
    case 'condfdr'  % include unconditioned FDR as additional column
        legends = cell(1,ntraits+1);
        ylabelstr = 'Conditional -log_{10}(FDR)';
        for i = 1:ntraits
            legends{i} = sprintf('%s | %s',traitname1,traitnames{i});
        end
        legends{end} = sprintf('%s',traitname1);
        
        fdrvec0( any(isnan(fdrmat),2) ) = NaN;
        logfdrmat = [logfdrmat, -log10(fdrvec0)];
        traitsequence = [ntraits+1, 1:ntraits];
        
    case 'conjfdr'
        legends = cell(1,ntraits);
        ylabelstr = 'Conjunctional -log_{10}(FDR)';
        for i = 1:ntraits
            legends{i} = sprintf('%s & %s',traitname1,traitnames{i});
        end
        traitsequence = 1:ntraits;
    otherwise
        error('')
end



%% PREPARE FIGURE -------------------------------------------------


% CORRECT the options stuff here

%handles = figure();
scrsz = get(0,'ScreenSize');
handles = figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)/2]);
figname = ['Manhattan ', sprintf('%s %s ',traitname1,traitnames{:}) ] ;
set(gca,'FontSize',20)
set(gcf,'Name',figname)
%set(gcf,'units','normalized','position',[0, 0, 1, 1]);
set(gca,'ColorOrder',colorlist)
if showgenes, ylim_extra = 3; else ylim_extra = 1; end
xlim([-100000,size(fdrmat,1)]), ylim([-0.1,max(logfdrmat(:))+ylim_extra ])
switch options.stattype
    case 'condfdr'
        colorlist(length(legends),:) = [0.5 0.5 0.5];  % black dots for unconditioned FDR
    case 'conjfdr'
end


%% DO MANHATTAN PLOT BY CHROMOSOME LIST

for chri = 1:length(chrnumlist)
    
    chr = chrnumlist{chri};
    ivec_chr = ismember(chrnumvec,chr);
    xtick = nan( length(chr), 1 );
    xticklabel = nan( length(chr), 1 );
    % Zoom on many or one chromosomes, alternating background color
    if length(chr)>1
        for k = chr
            ind = find(chrnumvec==k);
            if isempty(ind), continue; end;
            hold on;
            h=fill([ind(1) ind(end) ind(end) ind(1)],[0 0 300 300],[1 1 1] - graylevel * (1-mod(k,2)));
            set(h,'EdgeColor','none');
            if (~is_octave())
                set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            end
            xtick(k) = median(ind);
            xticklabel(k) = k;
        end
        set(gca,'XTick',xtick,'XTickLabel',xticklabel,'FontSize',fontsize_axes)
    else
        title(sprintf('Chromosome %d',chr),'FontSize',fontsize_axes);
        set(gca,'XTickLabel',[],'FontSize',fontsize_axes)
    end
    
    %  This is a dummy plot for legend purposes
    switch(options.stattype)
        case 'condfdr'
            legendsorder = [ntraits+1 , 1:ntraits];
        case 'conjfdr'
            legendsorder = 1:ntraits;
    end
    for i = legendsorder
        if ismember(i,traitsequence)
            hold on;
            plot(-1,-1,'o','MarkerFaceColor',colorlist(i,:),...
                'MarkerEdgeColor',colorlist(i,:),'MarkerSize',8);
        end
    end
    
    % Draw threshold line
    indvec = (1:length(chrnumvec))';
    hold on
    %plot([indvec(1), indvec(end)],-log10([options.fdrthresh, options.fdrthresh]),':k','LineWidth', 1)
    plot([-100000, indvec(end)],-log10([options.fdrthresh, options.fdrthresh]),':','color',[0,0,0],'LineWidth', 1.2)
    
    % Plot in three phases: (1) above threshold, (2) significant, (3) loci
    for phase=1:3
        for i = 1:size(logfdrmat,2)
            if ismember(i,traitsequence)
                if i<=ntraits
                    ivec = imat(:,i);
                    ivec2 = imat2(:,i);
                else % unconditioned FDR
                    ivec  = fdrvec0<=options.fdrthresh; % & isfinite(logfdrmat(:,i));
                    ivec2 = false(size(ivec));
                end
            end
            mrkface = colorlist(i,:);
            switch phase
                case 1   % do not pass threshold
                    selvec = ~ivec&ivec_chr;
                    mrkedge = colorlist(i,:);
                    mrksize = 4; lwid=1;
                case 2   % pass threshold, but not loci
                    selvec = ivec&~ivec2&ivec_chr;
                    mrkedge = colorlist(i,:);
                    mrksize = 8; lwid=1;
                case 3   % pass threshold and loci
                    selvec = ivec2&ivec_chr;
                    mrkedge = [0 0 0];
                    mrksize = 10; lwid=2;
            end
            hold on;
            [x, y] = filter_points_for_plotting(indvec(selvec), logfdrmat(selvec, i), options.manh_image_size);
            h=plot(x,y,'o','MarkerFaceColor',...
                mrkface,'MarkerEdgeColor',mrkedge,'LineWidth',lwid,'MarkerSize',mrksize);
            if (~is_octave())
                set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            end
        
            if(options.manh_redraw), drawnow, end
        end
    end
end

if length(chr)>1
    h=legend(legends(legendsorder),'Box','off','Location',LegendLocation);
    set(h,'FontSize',fontsize_legends);
end % Should automatically modify legends ith "trait1 & trait2"

switch options.stattype
    case 'condfdr'
        ylabel('-log_{10}(condFDR)','FontSize',fontsize_label);
    case 'conjfdr'
        ylabel('-log_{10}(conjFDR)','FontSize',fontsize_label);
end
xlabel('Chromosome','FontSize',fontsize_label);
%xtickangle(45)
set(gca,'LineWidth',1.2,'TickDir','out')
set(get(gca,'XAxis'),'TickLength',[0.007 0.007])
set(get(gca,'YAxis'),'TickLength',[0.004 0.004])

%% GENE NAMES -----------------------------------------
% One gene may have multiple loci according to LD
% Here we print one gene once, at one of the loci

if showgenes
    
    % Get snp-index (snpid) of unique genes (gnlist)
    idxloci = find( sum(imat2,2)>0 & ivec_chr & ~excludevec );    % snp at loci
    gnloci = genenamelist( idxloci );                % gene at loci
    if (is_octave())
        [gnlist,id,~] = unique(gnloci);              % list of unique genes
    else
        [gnlist,id,~] = unique(gnloci,'stable');     % list of unique genes
    end
    idxlist = idxloci(id);                          % SNP-idx of list of unique
    chrlist = chrnumvec(idxlist);                    % chromosome of SNP-idx
    
    % Print text boxes for gene names (textbox: [left,bottom,width,height])
    ymaxvec = max(logfdrmat,[],2); % max height of a snp for all traits
    ymaxvec(~ivec_chr) = 0;    
    textbox = NaN(length(gnlist),4);
    textboxhandle = NaN(length(gnlist),1);
    for i = 1:length(gnlist)
        xpos = idxlist(i);
        %ypos = ymaxvec(idxlist(i)) + ymargin;                          % fdr of single snp
        %ypos = max(ymaxvec) + ymargin;                                 % max fdr all
        %ypos = max( ymaxvec(chrnumvec==chrlist(i)) ) + ymargin;        % max fdr within chromosome
        
        % create box then reset y-position
        ypos = 0;
        h=text( xpos, ypos ,gnlist{i}, 'HorizontalAlignment','center',...
            'FontName','arial','FontSize',fontsize_genenames,'FontWeight','bold',...
            'FontAngle','italic','Color',[0.5 0.5 0.5]); %,'EdgeColor',[0,0,0]);
        textbox(i,:) = get(h,'Extent');
        x1 = max(floor(xpos-textbox(i,3)/2), 1);
        x2 = min(ceil(xpos+textbox(i,3)/2), size(imat2,1));
        textbox(i,2) = max( ymaxvec(x1:x2) ) + ymargin;
        set(h,'Position',[xpos,textbox(i,2)]);
        
        textboxhandle(i) = h;
    end
    if (options.manh_redraw), drawnow; end
    ysp = yspace*mean(textbox(:,4),1);    % white space
    
    % Rearrange text boxes position and color
    for i = 1:length(gnlist)   % plot from left to right
        
        % check if textbox collides with previous genenames
        left1   = textbox(i,1);
        right0  = textbox(1:i-1,1)+textbox(1:i-1,3);
        bottom1 = textbox(i,2);
        bottom0 = textbox(1:i-1,2);        
        
        % keep bumping up until it finds a space
        collide = left1<right0 & abs(bottom1-bottom0)< ysp;           
        while any(collide)
            %bottom1 = max(bottom0(left1<right0)) + ysp;
            %bottom1 = bottom0(end) + ysp;
            bottom1 = bottom1 + 0.25*ysp;
            textbox(i,2) = bottom1;
            collide = left1<right0 & abs(bottom1-bottom0)< ysp;
        end
    
        % Color according to trait with maximum conditional FDR
        % Special care should be taken to avoid coloring genes with a color of pruned trait.
        % For a given idxlist(i) we find which traits do not pass pruning
        % procedure, and exclude them from traitsequence. This works
        % for both conjfdr and condfdr.
        traitsequence_pruned = traitsequence(~ismember(traitsequence, find(~imat2(idxlist(i), :))));
        logfdrmat_pruned = logfdrmat(idxlist(i), traitsequence_pruned);
        [~, mi] = max(logfdrmat_pruned);
        mi = traitsequence_pruned(mi);
        xpos = textbox(i,1) + textbox(i,3)/2;
        ypos = textbox(i,2);
        set( textboxhandle(i), 'Position',[xpos,ypos],...
            'Color',colorlist(mi,:));      
    end
    
    % Rescale y-axis and gene fontsize
    if showgenes && ~isempty(textbox)
        yl_old = ylim();
        yl = [0, max(textbox(:,2))+ysp];
        ylim(yl)
        set( textboxhandle(:), 'FontSize', fontsize_genenames*yl_old(2)/yl(2) )
    end
end
