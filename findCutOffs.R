#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if ( length(args) < 1 ) {
  stop("I need a file with distances by many methods\n", call.=FALSE)
}
######### This script will cut a tree into groups

suppressWarnings(
    suppressPackageStartupMessages({
        library("tidyverse")
        library("cutpointr")
    })
)

## for plots:
plotFonts = theme(
    plot.title   = element_text(size = 33),
    axis.title.x = element_text(size = 26),
    axis.text.x  = element_text(size = 22),
    axis.title.y = element_text(size = 26),
    axis.text.y  = element_text(size = 22),
    )


fullCuts<-function(chidos,gachos,encabezado) {
    print(paste("calculating",chidos,"vs",gachos,encabezado,"cutoffs"),quote=F)

    tablita<-subset(tablota, Same == chidos | Same == gachos)
    medida<-get(encabezado)

    muchos=multi_cutpointr(tablita,
                           x=medidas,
                           class=Same,
                           direction = "<=",
                           pos_class=chidos,
                           neg_class=gachos,
                           method = maximize_metric,
                           metric = medida,
                           break_ties = median
                           )
    
    cuttable<-select(muchos,
                     optimal_cutpoint,
                     all_of(encabezado),
                     acc,
                     sensitivity,
                     specificity,
                     AUC,
                     predictor
                     )
    
    dashedtbl<-cuttable %>% mutate_if(is.character,
                                      str_replace_all,
                                      pattern = "\\.",
                                      replacement="-")
    
    roundedtbl<-dashedtbl %>% mutate_if(is.numeric, sprintf, fmt = "%.4f")
    
    outfile = paste(paste(chidos,gachos,sep="-"),
                    encabezado,"cutoffs",sep=".")
    if (!dir.exists(outdir)) {dir.create(outdir)}
    outfile<-paste(outdir,outfile,sep="/")
    write.table(roundedtbl,
                file=outfile,
                quote=FALSE,
                sep="\t",
                row.names=FALSE)

}

singleGraph<-function(ms,chidos,gachos) {
    print(paste("plotting",ms,"ROC:",chidos,"vs",gachos),quote=F)
    tablita<-subset(tablota, Same == chidos | Same == gachos)
    cuts=cutpointr(tablita,
                   x=!!ms,
                   class=Same,
                   direction = "<=",
                   pos_class=chidos,
                   neg_class=gachos,
                   method = maximize_metric,
                   metric = F1_score,
                   break_ties = median
                   )
    AUC=paste("AUC:",sprintf("%.3f",cuts$AUC))
    rocT<-gsub("\\.","-",ms)
    outfile<-paste(rocT,"pdf",sep=".")
    outfile<-paste(outdir,outfile,sep="/")
    roc<-plot_roc(cuts,display_cutpoint=F) + ggtitle(rocT) + geom_line(size=2) +
        theme_light() +
        annotate("rect", xmin = 0.2, xmax = 0.8, ymin = 0.4, ymax = 0.6,
                 fill="white",colour="black") +
        annotate("text", x = 0.5, y = 0.5, label = AUC, size = 14) +
        plotFonts
    ggsave(outfile)
}

outdir<-"CutoffTables"
print(paste("reading",args[1]),quote=F)
tablota<-read.table(args[1],sep="\t",head=T)
medidas<-grep("ANI|Signature|mash|dashing",colnames(tablota),value=T,perl=T)
metricas<-c("accuracy","F1_score","sum_sens_spec")
for ( metrica in metricas ) {
    fullCuts("Species","Genus",metrica)
    if( any(tablota$Same=="Family") ) {
        fullCuts("Species","Family",metrica)
        fullCuts("Genus","Family",metrica)
    }
}

for ( ms in medidas ) {
    singleGraph(ms,"Species","Genus")
}
