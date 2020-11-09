## Look at the properties for different intron types distinguished by their end sequences
## i.e. GTAG / ATAC, etc.

int.data <- read.table("Intron_ends_stats.csv", sep="\t", header=FALSE, stringsAsFactors=FALSE)

dim(int.data)
## [1] 4,368,399       6
## not too bad.
colnames(int.data) <- c('sp', 'db', 'gene', 'transcript', 'size', 'type')

dyn.load("strict_split/strict_split.so")
strict.split <- function(strings, split){
    .Call("strict_split", strings, split)
}

int.ts <- apply(int.data, 1, function(x){
    lst <- list('size'=as.numeric(strict.split(x['size'], ",")[[1]]),
                'type'=as.character(strict.split(x['type'], ",")[[1]]))
})

names(int.ts) <- int.data$gene

## let us harvest information about each species and each type of
## intron.

int.ts.t <- tapply( 1:length(int.ts), int.data$sp, function(i){
    x <- int.ts[i]
    ## get all sizes, all types and so on..
    sizes <- do.call(c, lapply(x, function(y) y$size ))
    types <- do.call(c, lapply(x, function(y) y$type ))
##
    ranks <- do.call(c, lapply(x, function(y){
        if(is.na( y$type[1] ))
            return(NA)
        1:length(y$type)
    }))
    ##
    row <- do.call(c, lapply(1:length(x), function(j){
        if(is.na( x[[j]]$type[1] ))
            return(NA)
        rep( i[j], length(x[[j]]$type )) }))
    list('i'=row, 'type'=types, 'size'=sizes, 'rank'=ranks)
})

## This looks good, so I do not know why I can't make dataframes of the
## above

sort(apply( t(sapply(int.ts.t, function(x){ sapply( x, length) })),
      1, sd ))

all.h <- hist( log2( 1 + unlist( lapply( int.ts, function(x){ x$size }))), breaks=20 )
    
type.size.h <- lapply( int.ts.t, function(x){
    if(length(x$size) != length(x$type))
        return(NA)
    b <- !is.na(x$size)
    sizes <- tapply( log2(x$size[b] + 1), x$type[b], function(x){ as.numeric(x) } )
    lapply( sizes, hist, plot=FALSE, breaks=all.h$breaks )
})






par(mfrow=c(2,2))
plot( type.size.h$homo_sapiens$GTAG )
plot( type.size.h$homo_sapiens$ATAC )

plot( type.size.h$danio_rerio$GTAG )
plot( type.size.h$danio_rerio$ATAC )

plot( type.size.h$takifugu_rubripes$GTAG )
plot( type.size.h$takifugu_rubripes$ATAC )

plot_hists <- function( x, n ){
    N <- sapply(x, function(y){ sum(y$counts) })
    o <- order(N, decreasing=TRUE)
    for(i in o[1:n]){
        plot(x[[i]], main=names(x)[i] )
    }
}

par(mfrow=c(3, 4))
plot_hists( type.size.h$homo_sapiens, n=12 )
plot_hists( type.size.h$mus_musculus, n=12 )
plot_hists( type.size.h$danio_rerio, n=12 )

par(mfrow=c(3, 4))
for(sp in names(type.size.h)){
    plot_hists( type.size.h[[sp]], n=12 )
    inpt <- readline(paste(sp, " "))
}


## let us do this for
sel.sp <- c('GTAG', 'GCAG', 'ATAC')

for(sp in names(type.size.h)){
    if( all( sel.sp %in% names( type.size.h[[sp]] ))){
        density <- sapply( type.size.h[[sp]][sel.sp], function(x){ x$density })
        plot( all.h$mids, density[,1], ylim=range(density), type='n', main=sp )
        for(i in 1:length(sel.sp)){
            lines( all.h$mids, density[,i], col=i, lwd=2 )
        }
        legend('topright', sel.sp, col=1:length(sel.sp), lty=1, lwd=2)
        inpt <- readline(paste(sp, " "))
    }
}
