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
gene_structure_summary <- function(utr5_lengths, cds_lengths, utr3_lengths) {
    # equalize lengths
    max_length <- max(length(utr5_lengths), length(cds_lengths), length(utr3_lengths))
    utr5_lengths <- c(utr5_lengths, rep(NA, max_length - length(utr5_lengths)))
    cds_lengths  <- c(cds_lengths,  rep(NA, max_length - length(cds_lengths)))
    utr3_lengths <- c(utr3_lengths, rep(NA, max_length - length(utr3_lengths)))

    mat <- cbind(utr5_lengths, cds_lengths, utr3_lengths)

    summary_dat <- data.frame(
        min=apply(mat, 2, function (x) { min(x, na.rm=TRUE) }),
        median=apply(mat, 2, function (x) { median(x, na.rm=TRUE) }),
        mean=apply(mat, 2, function (x) { mean(x, na.rm=TRUE) }),
        max=apply(mat, 2, function (x) { max(x, na.rm=TRUE) })
    )
    rownames(summary_dat) <- c("5'UTR", "CDS", "3'UTR")
    
    return(summary_dat)
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

