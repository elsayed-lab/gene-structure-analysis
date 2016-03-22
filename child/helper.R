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
