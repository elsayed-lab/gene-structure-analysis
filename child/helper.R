#
# plot_gene_structure
#
# based on:
# http://stackoverflow.com/questions/9100841/rectangle-bar-graph-filled-with-color-and-distance-using-r-base-r-or-ggplot2-or
#
plot_gene_structure <- function(utr5, cds, utr3, utr_color='#6a839c', cds_color='#4b739c') {
    # determine boundaries
    xboundaries <- c(0, utr5, utr5 + cds, utr5 + cds + utr3)
    xticks <- c((utr5/2),
                (cds /2) + utr5,
                (utr3/2) + utr5 + cds)

    gene_structure_labels <- c(sprintf("5'UTR (%dnt)", round(utr5)),
                               sprintf("CDS (%dnt)", round(cds)),
                               sprintf("3'UTR (%dnt)", round(utr3)))
    gene_structure_colors <- c(utr_color, cds_color, utr_color)

    # create plot
    barplot(as.matrix(diff(xboundaries)), horiz=T, col=gene_structure_colors, axes=F)
    axis(1, labels=gene_structure_labels, at=xticks)
}

#
# gene_structure_summary
#
gene_structure_summary <- function(lengths) {
    data.frame(
        min=apply(lengths, 2, function (x) { min(x, na.rm=TRUE) }),
        median=apply(lengths, 2, function (x) { median(x, na.rm=TRUE) }),
        mean=apply(lengths, 2, function (x) { mean(x, na.rm=TRUE) }),
        max=apply(lengths, 2, function (x) { max(x, na.rm=TRUE) })
    )
}

#
# get_intergenic_lengths
#
get_intergenic_lengths <- function(gr, min_read_support=1) {
    # Vector to keep track of intergenic lengths
    result <- data.frame()

    # Make sure the GenomicRanges instance is in sorted order
    gr <- sortSeqlevels(gr)
    gr <- sort(gr, ignore.strand=TRUE)

    # UTR GFF feature annotations
    utr_types <- c('five_prime_UTR', 'three_prime_UTR')

    # Iterate over chromosome
    for (ch in levels(seqnames(gr))) {
        # Get GRange subsets for the chromosome
        gr_ch <- gr[seqnames(gr) == ch]

        # Iterate over features on chromosome
        # we can skip the last feature since it cannot fall between two genes
        # on the same strand
        for (i in 1:(length(gr_ch) - 1)) {
            left  <- gr_ch[i]
            right <- gr_ch[i + 1]

            # Check for adjacent UTRs
            # TODO: add max distance to exclude inter-PTU regions
            if (left$type %in% utr_types && right$type %in% utr_types) {
                # Check to make sure lefts have minimum required read
                # support
                if (left$score < min_read_support || right$score < min_read_support) {
                    next
                }
                
                # Skip cases where type is same (can only occur at switch
                # sites).
                if (left$type == right$type) {
                    next
                }
                intergenic_length <- start(right) - end(left)
                result <- rbind(result, c(left$ID, right$ID, intergenic_length))
                message(sprintf("%s (i=%d): %d", ch, i, intergenic_length))
            }
        }
    }
    colnames(result) <- c("left", "right", "length")
    result$length <- as.numeric(result$length)

    return(result)
}

#
# get_intercds_sequences
#
# Finds the inter-CDS sequence for each pair of adjacent CDS's on the same
# strand
#
get_intercds_sequences <- function(gr, fasta) {
    # Vector to keep track of inter-CDS sequences
    result <- c()

    # Make sure the GenomicRanges instance is in sorted order
    gr <- sortSeqlevels(gr)
    gr <- sort(gr)

    # Iterate over chromosome
    for (ch in levels(seqnames(gr))) {
        # Get GRange subsets for the chromosome
        gr_ch <- gr[seqnames(gr) == ch]

        gr_cds <- gr_ch[gr_ch$type == 'CDS']

        # Iterate over CDS's on chromosome
        # we can skip the last feature since it cannot fall between two genes
        # on the same strand
        for (i in 1:(length(gr_cds) - 1)) {
            cds <- gr_cds[i]
            next_cds <- gr_cds[i + 1]

            if (strand(cds) == strand(next_cds)) {
                result <- append(result, subseq(fasta[ch], end(cds), start(next_cds)))
            }
        }
    }

    return(result)
}


#
# detect_polypyrimidine_tracts
#
detect_polypyrimidine_tracts <- function() {

}

#' 
#' create_utr_comparison_df
#' 
create_utr_comparison_df <- function(sites, s1_name, s2_name, min_read_support) {
    # Get sites for specified developmental stages
    s1 <- sites[[s1_name]] %>% 
        filter(type=='primary') %>%
        arrange(gene)
    s2 <- sites[[s2_name]] %>% 
        filter(type=='primary') %>%
        arrange(gene)

    # limit to genes for which primary acceptor site was detected in both stages
    common_genes <- intersect(s1$gene, s2$gene)

    s1 <- s1 %>% filter(gene %in% common_genes)
    s2 <- s2 %>% filter(gene %in% common_genes)

    # for each gene, get maximum coverage for alternative sites; used to compute
    # primary/secondary ratio
    s1_alt_coverage <- sites[[s1_name]] %>% 
        filter(type=='alternative') %>% 
        group_by(gene) %>% 
        top_n(1, num_reads) %>% 
        slice(1) %>%
        select(gene, alt_coverage=num_reads) %>%
        ungroup()

    # default to alt coverage of "1" (neccessary to compute ratios)
    s1_no_alt_sites <- s1$gene[!s1$gene %in% s1_alt_coverage$gene]
    s1_alt_coverage <- rbind(s1_alt_coverage,
                             data.frame(gene=s1_no_alt_sites, alt_coverage=1))

    s1 <- merge(s1, s1_alt_coverage, by='gene') %>%
        mutate(ptos=num_reads / alt_coverage)

    s2_alt_coverage <- sites[[s2_name]] %>% 
        filter(type=='alternative') %>% 
        group_by(gene) %>% 
        top_n(1, num_reads) %>% 
        slice(1) %>%
        select(gene, alt_coverage=num_reads) %>%
        ungroup()

    # default to alt coverage of "1" (neccessary to compute ratios)
    s2_no_alt_sites <- s2$gene[!s2$gene %in% s2_alt_coverage$gene]
    s2_alt_coverage <- rbind(s2_alt_coverage,
                             data.frame(gene=s2_no_alt_sites, alt_coverage=1))

    s2 <- merge(s2, s2_alt_coverage, by='gene') %>%
        mutate(ptos=num_reads / alt_coverage)

    # combine into a single dataframe (both data frames have been previously sorted
    # by gene id to ensure correct order)
    dat <- data_frame(
        gene=s1$gene,
        s1_len=s1$utr_length,
        s2_len=s2$utr_length,
        s1_num_reads=s1$num_reads,
        s2_num_reads=s2$num_reads,
        s1_ptos=s1$ptos,
        s2_ptos=s2$ptos
    )

    # ignore UTRs of the same length (uninteresting) and add difference column
    dat <- dat %>% 
        filter(s1_len != s2_len) %>% 
        mutate(len_diff=abs(s1_len - s2_len))

    # use min read support as a measure of our confidence
    dat <- dat %>% 
        mutate(min_num_reads=pmin(s1_num_reads, s2_num_reads))
    #    mutate(log_min_num_reads=log(min_num_reads))

    # only show sites with > 3 reads in both stages
    dat %>% filter(min_num_reads >= min_read_support)
}

plot_diff_utrs <- function(dat, feature_name) {
    # plot UTR length differences across developmental stages
    ggplot(dat, aes(s1_len, s2_len, color=min_num_reads, size=min_num_reads)) + 
        geom_abline(slope=1, intercept=0, color='#CCCCCC', lwd=0.5) +
        geom_abline(slope=1, intercept=300, color='#666666', lwd=0.5) +
        geom_abline(slope=1, intercept=-300, color='#666666', lwd=0.5) +
        geom_point() +
        xlab(s1_name) + ylab(s2_name) +
        scale_color_gradient2(low="black", mid="blue", high="red") + 
        scale_size_continuous() +
        scale_x_continuous(limits=c(0,2000), expand=c(0.01, 0.01)) +
        scale_y_continuous(limits=c(0,2000), expand=c(0.01, 0.01)) + 
        #geom_text(data=subset(dat, min_num_reads >= 100 & len_diff > 300),
        #          aes(s1_len, s2_len, label=gene),
        #          color='#000000', size=4, angle=45, hjust=0, vjust=0) +
        ggtitle(sprintf("%s site usage: %s vs. %s (All genes)", feature_name, s1_name, s2_name)) +
        theme(text=element_text(size=12, family='DejaVu Sans'),
            plot.title=element_text(size=rel(1)))
}

#'
#' plot_alt_site_distance_hist
#'
plot_alt_site_distance_hist <- function(sites, gene_strands,
                                        upstream_color='red',
                                        downstream_color='blue',
                                        stroke_color='#333333', xlabel='',
                                        main='') {
    # get primary and alternative sites
    primary_sites <- sites %>% 
        filter(type=='primary') %>%
        select(gene, chr_coordinate) %>%
        rename(loc_primary=chr_coordinate) %>%
        arrange(gene)

    secondary_sites <- sites %>% 
        filter(type=='alternative') %>%
        select(gene, chr_coordinate) %>%
        rename(loc_secondary=chr_coordinate) %>%
        arrange(gene)

    # exclude genes with no detected alternative sites
    primary_sites <- primary_sites %>% filter(gene %in% secondary_sites$gene)

    # combine into a single dataframe
    dat <- merge(primary_sites, secondary_sites, by='gene', all=TRUE)

    # add strand information
    dat <- merge(dat, gene_strands)

    # determine distance and direction of secondary sites relative to primary ones
    dat <- dat %>% mutate(dist=ifelse(strand == '+', 
                                    loc_secondary - loc_primary,
                                    loc_primary - loc_secondary)) %>%
                mutate(direction=factor(sign(dist), 
                                        labels=c('upstream', 'downstream')))

    # limit to sites in range [-2500, 2500]
    dat <- dat %>% filter(abs(dist) <= 2500)

    ggplot(dat, aes(dist, fill=direction)) +
        geom_histogram(breaks=seq(-2500, 2500, by=50), color=stroke_color) +
        scale_fill_manual(values=c(upstream_color, downstream_color)) +
        ggtitle(main) + 
        theme(plot.title=element_text(hjust=0)) + 
        labs(x=xlabel)
}

