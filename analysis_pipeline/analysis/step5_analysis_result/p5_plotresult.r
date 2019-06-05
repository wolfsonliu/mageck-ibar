args <-  commandArgs(trailingOnly=TRUE)

library(ggplot2)

## args <- c()
## args[1] <- '/gpfs/share/home/1501111485/Project/Little/20190102/data/ibar/A/iBAR_A'
## args[2] <- '/gpfs/share/home/1501111485/Project/Little/20190102/process/p5_plotresult/fig/iBAR_A'

filenames <- c(
    'gene.high'=paste(args[1], 'gene.high.txt', sep='.'),
    'gene.low'=paste(args[1], 'gene.low.txt', sep='.')
)

data <- list()

for (x in names(filenames)) {
    data[[x]] <- read.table(
        filenames[[x]], header=TRUE, stringsAsFactors=FALSE
    )
    data[[x]]$type <- 'Gene'
    data[[x]]$type[grepl('neg', data[[x]]$group_id)] <- 'Non-targeting'
    data[[x]]$type[
                  order(data[[x]]$lo_value, decreasing=FALSE)[1:10]
              ] <- data[[x]]$group_id[
                                 order(data[[x]]$lo_value, decreasing=FALSE)[1:10]
                             ]
    data[[x]]$type <- factor(
        data[[x]]$type,
        levels=c(
            'Gene', 'Non-targeting',
            data[[x]]$group_id[order(data[[x]]$lo_value, decreasing=FALSE)[1:10]]
        )
    )
    plotcolors <- c(
        'black', 'lightgray',
        '#e31a1c', '#ff7f00', '#1f78b4', '#33a02c', '#6a3d9a',
        '#fb9a99', '#fdbf6f', '#a6cee3', '#b2df8a', '#cab2d6'
    )
    names(plotcolors) <- c(
        'Non-targeting', 'Gene',
        data[[x]]$group_id[order(data[[x]]$lo_value, decreasing=FALSE)[1:10]]
    )
    data[[x]]$size <- 1.1
    data[[x]]$size[data[[x]]$type == 'Gene'] <- 0.9
    data[[x]]$size[data[[x]]$type == 'Non-targeting'] <- 0.8
    data[[x]]$x <- sample(
        seq(round(dim(data[[x]])[1] / 5)), dim(data[[x]])[1], replace=TRUE
    )

    p <- ggplot(
        data[[x]][order(data[[x]]$type),],
        aes(x=x, y=-log10(lo_value), color=type, size=size)
    ) + geom_point(
            alpha=0.8
        ) + scale_color_manual(
                values=plotcolors
            ) + theme(
                    panel.background=element_rect(fill='transparent', color='black'),
                    panel.grid.major.x=element_blank(),
                    panel.grid.minor.x=element_blank(),
                    panel.grid.major.y=element_line(
                        color="#d9d9d9", linetype='dashed'
                    ),
                    panel.grid.minor.y=element_blank(),
                    axis.text.x=element_blank(),
                    axis.ticks.x=element_blank()
                ) + xlab(
                        ''
                    ) + ylab(
                            expression(paste('RRA Score (', -log[10], beta, ')'))
                        ) + guides(size = "none")
    ggsave(
        paste(args[2], 'screenscatter', x, 'pdf', sep='.'), p
    )
}
