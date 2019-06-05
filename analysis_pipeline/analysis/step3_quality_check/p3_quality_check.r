library(ggplot2)
library(dplyr)
library(tidyr)

dirs <- list()
dirs[['base']] <- '/gpfs/share/home/1501111485/Project/Little/20190102'
dirs[['data']] <- file.path(dirs[['base']], 'data')
dirs[['count']] <- file.path(dirs[['data']], 'count')
dirs[['analysis']] <- file.path(dirs[['base']], 'process', 'p3_quality_check')
dirs[['fig']] <- file.path(dirs[['analysis']], 'fig')

## A_R1_S2_L003_R1_001.fastq.gz: 86047592
## A_R2_S3_L003_R1_001.fastq.gz: 73534950
## Ctr_R1_S4_L003_R1_001.fastq.gz: 40128871
## Ctr_R2_S5_L003_R1_001.fastq.gz: 60722364
## D1_R1_S6_L003_R1_001.fastq.gz: 19882938
## D1_R2_S7_L003_R1_001.fastq.gz: 29736320
## D4_R2_S6_L004_R1_001.fastq.gz: 52366204
## M_R1_S3_L003_R1_001.fastq.gz: 72566074
## M_R2_S4_L003_R1_001.fastq.gz: 118280513
## U1_R1_S3_L004_R1_001.fastq.gz: 46859352
## U1_R2_S4_L004_R1_001.fastq.gz: 31917456
## U4_R1_S1_L003_R1_001.fastq.gz: 15993197
## U4_R2_S2_L003_R1_001.fastq.gz: 16847653

seqcounts <- c('A_R1'=86047592, 'A_R2'=73534950, 'Ctr_R1'=40128871,
               'Ctr_R2'=60722364, 'D1_R1'=19882938, 'D1_R2'=29736320,
               'D4_R2'=52366204, 'M_R1'=72566074, 'M_R2'=118280513,
               'U1_R1'=46859352, 'U1_R2'=31917456, 'U4_R1'=15993197,
               'U4_R2'=16847653)

## A_R1.rawcount: 44708092
## A_R2.rawcount: 38598102
## Ctr_R1.rawcount: 21453801
## Ctr_R2.rawcount: 31764090
## D1_R1.rawcount: 10505754
## D1_R2.rawcount: 15719200
## D4_R2.rawcount: 25713463
## M_R1.rawcount: 34913082
## M_R2.rawcount: 58952483
## U1_R1.rawcount: 24902619
## U1_R2.rawcount: 17097036
## U4_R1.rawcount: 8666941
## U4_R2.rawcount: 9024631

rawcounts <- c('A_R1'=44708092, 'A_R2'=38598102, 'Ctr_R1'=21453801,
               'Ctr_R2'=31764090, 'D1_R1'=10505754, 'D1_R2'=15719200,
               'D4_R2'=25713463, 'M_R1'=34913082, 'M_R2'=58952483,
               'U1_R1'=24902619, 'U1_R2'=17097036, 'U4_R1'=8666941,
               'U4_R2'=9024631)

libs <- list()
for (x in c('A', 'D1', 'M', 'U1', 'U4')) {
    libs[[x]] <- read.csv(
        file.path(dirs[['count']], paste(x, 'count', 'csv', sep='.')),
        header=TRUE, stringsAsFactors=FALSE
    )
}

libcounts <- c(
    'A_R1'=sum(libs[['A']][['A_R1']]),
    'A_R2'=sum(libs[['A']][['A_R2']]),
    'Ctr_R1'=sum(libs[['A']][['Ctr_R1']]),
    'Ctr_R2'=sum(libs[['A']][['Ctr_R2']]),
    'D1_R1'=sum(libs[['D1']][['D1_R1']]),
    'D1_R2'=sum(libs[['D1']][['D1_R2']]),
    'D4_R2'=sum(libs[['D4']][['D4_R2']]),
    'M_R1'=sum(libs[['M']][['M_R1']]),
    'M_R2'=sum(libs[['M']][['M_R2']]),
    'U1_R1'=sum(libs[['U1']][['U1_R1']]),
    'U1_R2'=sum(libs[['U1']][['U1_R2']]),
    'U4_R1'=sum(libs[['U4']][['U4_R1']]),
    'U4_R2'=sum(libs[['U4']][['U4_R2']])
)

countdata <- as.data.frame(cbind(seqcounts, rawcounts, libcounts))
colnames(countdata) <- c('Read.Counts', 'Raw.Counts', 'Library.Counts')
countdata$Exp. <- rownames(countdata)

countdata.rate <- countdata
countdata.rate$Read.Counts <- 1
countdata.rate$Raw.Counts <- countdata$Raw.Counts / countdata$Read.Counts
countdata.rate$Library.Counts.1 <- countdata$Library.Counts / countdata$Read.Counts
countdata.rate$Library.Counts.2 <- countdata$Library.Counts / countdata$Raw.Counts
countdata.rate$Library.Counts <- NULL

pdata <- gather(countdata, key='Type', value='Counts', -Exp.)
pdata$Type <- factor(
    pdata$Type, levels=c('Read.Counts', 'Raw.Counts', 'Library.Counts'),
    labels=c('Read Counts', 'Raw Counts', 'Library Counts')
)

pdata.rate <- pdata[pdata$Type != 'Read Counts',]
pdata.rate <- left_join(
    pdata.rate, pdata[pdata$Type == 'Read Counts', c('Exp.', 'Counts')],
    by=c('Exp.'), suffix=c('', '.Read')
)
pdata.rate <- left_join(
    pdata.rate, pdata[pdata$Type == 'Raw Counts', c('Exp.', 'Counts')],
    by=c('Exp.'), suffix=c('', '.Raw')
)
pdata.rate$Rate.Read <- pdata.rate$Counts / pdata.rate$Counts.Read
pdata.rate$Rate.Raw <- pdata.rate$Counts / pdata.rate$Counts.Raw
pdata.rate$Rate.Read.t <- paste0(round(pdata.rate$Rate.Read * 100), '%')
pdata.rate$Rate.Raw.t <- paste0(round(pdata.rate$Rate.Raw * 100), '%')
pdata.rate <- pdata.rate[pdata.rate$Counts != 0,]
pdata.rate$Counts.Read <- NULL
pdata.rate$Counts.Raw <- NULL

p <- ggplot(
    pdata, aes(x=Exp., y=Counts, group=Type, fill=Type)
) + geom_bar(
        stat='identity', position='dodge'
    ) + geom_text(
            data=pdata.rate[pdata.rate$Type == 'Raw Counts',],
            mapping=aes(x=Exp., y=Counts, label=paste0(Rate.Read.t, '\n')),
            nudge_x=0, nudge_y=0.
        ) + geom_text(
                data=pdata.rate[pdata.rate$Type == 'Library Counts',],
                mapping=aes(
                    x=Exp., y=Counts,
                    label=paste(Rate.Read.t, Rate.Raw.t, '\n', sep='\n')
                ),
                nudge_x=0.3, nudge_y=0.8
            ) + theme(
                    panel.background=element_rect(fill='transparent', color='black'),
                    panel.grid.major.x=element_blank(),
                    panel.grid.minor.x=element_blank(),
                    panel.grid.major.y=element_line(
                        color="#d9d9d9", linetype='dashed'
                    ),
                    panel.grid.minor.y=element_blank()
                ) + xlab('Experiment') + ylab('Counts')

ggsave(file.path(dirs[['fig']], 'count.pdf'), p, width=40, height=10, units='cm')

####################
## Plot
####################

libs.norm <- libs

for (x in names(libs)) {
    for (y in colnames(libs[[x]])[-c(1,2,3)]) {
        libs.norm[[x]][[y]] <- libs[[x]][[y]] / sum(libs[[x]][[y]]) * 10^6
    }
}

## Boxplot
pdata <- unique(Reduce(
    rbind,
    lapply(
        libs.norm,
        function(x) {
            gather(
                x, key='Exp.', value='RPM', -c('gene', 'guide', 'barcode')
            )
        }
    )
))

pdata$Exp. <- factor(
    pdata$Exp.,
    levels=c('Ctr_R1', 'Ctr_R2', 'A_R1', 'A_R2',  'M_R1', 'M_R2',
             'D1_R1', 'D1_R2', 'U1_R1', 'U1_R2', 'U4_R1', 'U4_R2')
)

p <- ggplot(
    pdata, aes(x=Exp., y=RPM + 1,)
) + geom_boxplot(
    ) + theme(
            panel.background=element_rect(fill='transparent', color='black'),
            panel.grid.major.x=element_blank(),
            panel.grid.minor.x=element_blank(),
            panel.grid.major.y=element_line(
                color="#d9d9d9", linetype='dashed'
            ),
            panel.grid.minor.y=element_blank()
        ) + scale_y_log10() + xlab('Experiment') + ylab('Reads Per Million + 1')

ggsave(file.path(dirs[['fig']], 'boxplot.pdf'), p, width=20, height=10, units='cm')

## Scatterplot
pdata <- unique(Reduce(
    rbind,
    lapply(
        libs.norm,
        function(x) {
            gather(
                x, key='Exp.', value='RPM', -c('gene', 'guide', 'barcode')
            )
        }
    )
))
pdata$Group <- unlist(lapply(strsplit(pdata$Exp., '_'), '[', 1))
pdata$Rep. <- unlist(lapply(strsplit(pdata$Exp., '_'), '[', 2))
pdata$Exp. <- NULL
pdata <- spread(pdata, key='Rep.', value='RPM')

pdata$type <- pdata$Group
pdata$type[grepl('negative', pdata$gene)] <- 'Non-targeting'
pdata$type <- factor(
    pdata$type,
    levels=c('Non-targeting', 'Ctr', 'A', 'D1', 'M', 'U1', 'U4')
)

pdata$Group <- factor(
    pdata$Group,
    levels=c('Ctr', 'A', 'D1', 'M', 'U1', 'U4')
)

p <- ggplot(
    pdata, aes(x=R1, y=R2, color=type)
) + geom_point(
    ) + scale_color_manual(
            values=c(
                'Non-targeting'='black',
                'Ctr'='#a65628', 'A'='#e41a1c', 'D1'='#377eb8', 'M'='#4daf4a',
                'U1'='#984ea3', 'U4'='#ff7f00'
            )
        ) + facet_wrap(
                .~Group
            ) + theme(
                    panel.background=element_rect(fill='transparent', color='black'),
                    panel.grid.major.x=element_line(
                        color="#d9d9d9", linetype='dashed'
                    ),
                    panel.grid.minor.x=element_blank(),
                    panel.grid.major.y=element_line(
                        color="#d9d9d9", linetype='dashed'
                    ),
                    panel.grid.minor.y=element_blank()
                ) + scale_x_log10(
                    ) + scale_y_log10(
                        ) + xlab(
                                'Replicate 1 (RPM + 1)'
                            ) + ylab(
                                    'Replicate 2 (RPM + 1)'
                                ) + coord_fixed(ratio=1)

ggsave(file.path(dirs[['fig']], 'scatterplot.pdf'), p, width=40, height=10, units='cm')

## Scatterplot separated

for (x in unique(pdata$Group)) {
    p <- ggplot(
        pdata[pdata$Group == x,], aes(x=R1, y=R2, color=type)
    ) + geom_point(
        ) + scale_color_manual(
                values=c(
                    'Non-targeting'='black',
                    'Ctr'='#a65628', 'A'='#e41a1c', 'D1'='#377eb8', 'M'='#4daf4a',
                    'U1'='#984ea3', 'U4'='#ff7f00'
                )
            ) + facet_grid(
                    .~Group
                ) + theme(
                        panel.background=element_rect(fill='transparent', color='black'),
                        panel.grid.major.x=element_line(
                            color="#d9d9d9", linetype='dashed'
                        ),
                        panel.grid.minor.x=element_blank(),
                        panel.grid.major.y=element_line(
                            color="#d9d9d9", linetype='dashed'
                        ),
                        panel.grid.minor.y=element_blank()
                    ) + scale_x_log10(
                        ) + scale_y_log10(
                            ) + xlab(
                                    'Replicate 1 (RPM + 1)'
                                ) + ylab(
                                        'Replicate 2 (RPM + 1)'
                                    ) + coord_fixed(ratio=1)
    ggsave(
        file.path(dirs[['fig']], paste('scatterplot', x, 'pdf', sep='.')), p
    )
}
