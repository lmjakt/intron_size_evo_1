## This redoes many of the analyses done in ../R_Mosjoen/
## but with a larger number of species. Insetead of 20 vertebrates most of which are fish
## we have 172, most of which are mammals..


source("../R_commmon_functions/functions.R")
source("~/R/general_functions.R")
require(RMySQL)

gene.stats <- read.table("../family_members/orthologue_transcripts/exon_intron_stats.csv", sep="\t", stringsAsFactors=FALSE)
colnames(gene.stats) <- c('sp', 'db', 'family', 'gene', 'transcript', 'chr', 'strand', 'pos', 'ex.s', 'in.s')

dim(gene.stats)
## [1] 1063202      10
## about a million rows

t.lengths <- sapply( gene.stats[,'ex.s'], function(x){
    sum( as.numeric( strsplit( x, ',' )[[1]] ) ) })

sum(t.lengths > 10000) / length(t.lengths)
## [1] 0.006477603

## how many families have a maximum length
## could probably have done this using the t.lengths structure.. 
fam.tlengths <- tapply( 1:nrow(gene.stats), gene.stats[,'family'], function(i){
    sapply( gene.stats[i,'ex.s'], function(x){
        sum( as.numeric( strsplit( x, ',' )[[1]] )) })
})

fam.max.tlengths <- sapply( fam.tlengths, max )
sum( fam.max.tlengths > 10000 ) ## 506
sum( fam.max.tlengths > 15000 ) / length( fam.max.tlengths) ## [1] 0.01193981

## the above are to avoid very long sequenc alignments which take up
## too much memory.

dbs <- strsplit(readLines( "../family_members/vertebrate_family_members_1_130.txt", n=1 ), "\t")[[1]]
names(dbs) <- sub("_core_98_\\d+", "", dbs )

length(dbs)
## 172

genome.stats <- lapply(dbs, get.genome.stats)
genome.sizes <- sapply(genome.stats, function(x){ x[ x$statistic == 'ref_length', 'value' ] })
write.table(genome.sizes, 'genome_sizes.txt', row.names=TRUE, quote=FALSE)
fam.sizes <- tapply( 1:nrow(gene.stats), gene.stats$family, length )

sp.gene.n <- tapply( 1:nrow(gene.stats), gene.stats$sp, length )

par('mar'=c(8.1, 4.1, 4.1, 2.1))
barplot( sp.gene.n, names.arg=sub("_[^_]+$", "", names(sp.gene.n)), las=2)

exon.s <- lapply( strsplit( gene.stats$ex.s, "," ), as.numeric )
intron.s <- lapply( strsplit( gene.stats$in.s, ","), as.numeric )

intron.d.all <- hist( log2( unlist( intron.s )), breaks=40 )

intron.d <- tapply( intron.s, gene.stats$sp, function(x){
    hist( log2( unlist( x ) ), breaks=intron.d.all$breaks, plot=FALSE)
})

intron.d1 <- tapply( intron.s, gene.stats$sp, function(x){
    hist( log2(unlist( sapply( x, function(y){ y[1] } ))), breaks=intron.d.all$breaks, plot=FALSE )
})

## an ordered list of species names
sp.n <- names( sort(genome.sizes) )

##pdf("intron_size_distributions_raw.pdf", width=16, height=12)
par(mfrow=c(4,5))
par(oma=c(0,0,0,0))
par(mar=c(5.1,4.1,4.1,2.1))
for( b in seq(1, length(sp.n), 20 )){
    for(name in sp.n[b + 0:19]){
        plot( intron.d[[name]], main=sub("_", " ", name), xlab='log2 intron length' )
    }
    inpt <- readline("next")
}
##dev.off()

## some of these should clearly be removed from the analysis, due to bad annotation
##
### try a violin plot instead..
## use density
intron.dd <- tapply( intron.s, gene.stats$sp, function(x){
    density( log2( unlist(x) ), na.rm=TRUE)
})

intron.dd.1 <- tapply( intron.s, gene.stats$sp, function(x){
    density( log2(unlist( sapply( x, function(y){ y[1] } ))), na.rm=TRUE )
})

intron.dd.2 <- tapply( intron.s, gene.stats$sp, function(x){
    density( log2(unlist( sapply( x, function(y){ y[2] } ))), na.rm=TRUE )
})

intron.dd.3 <- tapply( intron.s, gene.stats$sp, function(x){
    density( log2(unlist( sapply( x, function(y){ y[3] } ))), na.rm=TRUE )
})

par(mfrow=c(1,1))
par('mar'=c(7.1, 4.1, 4.1, 5.1))
plot.new()
plot.window(xlim=c(0.5, 0.5+length(intron.dd)),
            ylim=range( c(unlist(sapply(intron.dd, function(x){ x$x })) )),
            yaxs='i', xaxs='i' )
tmp.gm <- 5 / max(genome.sizes)
tmp.d1 <- unlist(sapply( intron.dd, function(x){ x$y }))
tmp.d2 <- unlist(sapply( intron.dd.1, function(x){ x$y }))
tmp.m <- 0.5 / max(c(tmp.d1, tmp.d2))
for(i in 1:length(sp.n)){
    name <- sp.n[i]
    x.d1 <- intron.dd[[name]]$y * tmp.m
    x.d2 <- intron.dd.1[[name]]$y * tmp.m
    x1 <- c(i-x.d1, i+rev(x.d1))
    x2 <- c(i-x.d2, i+rev(x.d2))
    y1 <- c(intron.dd[[name]]$x, rev(intron.dd[[name]]$x))
    y2 <- c(intron.dd.1[[name]]$x, rev(intron.dd.1[[name]]$x))
    rect(i-0.2, 0, i+0.2, genome.sizes[name] * tmp.gm, col=rgb(0.4, 0.4, 0.4, 0.4), border=NA)
    polygon(x1,y1, col=rgb(0, 0.1, 0.8, 0.3), border=NA)
##    polygon(x2,y2, col=rgb(0.8, 0.1, 0, 0.3), border=NA)
}
axis(side=1, at=1:length(intron.dd), sub("_[^_]+$", "", sp.n), las=2)
axis(side=2, at=seq(0, 20, 4))
tmp <- seq(0, 3e9, 1e9)
axis(side=4, at=tmp*tmp.gm, labels=tmp, line=0, las=2)
mtext("log 2 intron length", 2, line=3)

## include first introns separately

par(mfrow=c(1,1))
par('mar'=c(7.1, 4.1, 4.1, 5.1))
plot.new()
plot.window(xlim=c(0.5, 0.5+length(intron.dd)),
            ylim=range( c(unlist(sapply(intron.dd, function(x){ x$x })) )),
            yaxs='i', xaxs='i' )
tmp.gm <- 5 / max(genome.sizes)
tmp.d1 <- unlist(sapply( intron.dd, function(x){ x$y }))
tmp.d2 <- unlist(sapply( intron.dd.1, function(x){ x$y }))
tmp.m <- 0.5 / max(c(tmp.d1, tmp.d2))
for(i in 1:length(sp.n)){
    name <- sp.n[i]
    x.d1 <- intron.dd[[name]]$y * tmp.m
    x.d2 <- intron.dd.1[[name]]$y * tmp.m
    x1 <- c(i-x.d1, i+rev(x.d1))
    x2 <- c(i-x.d2, i+rev(x.d2))
    y1 <- c(intron.dd[[name]]$x, rev(intron.dd[[name]]$x))
    y2 <- c(intron.dd.1[[name]]$x, rev(intron.dd.1[[name]]$x))
    rect(i-0.2, 0, i+0.2, genome.sizes[name] * tmp.gm, col=rgb(0.4, 0.4, 0.4, 0.4), border=NA)
    polygon(x1,y1, col=rgb(0, 0.1, 0.8, 0.3), border=NA)
    polygon(x2,y2, col=rgb(0.8, 0.1, 0, 0.3), border=NA)
}
axis(side=1, at=1:length(intron.dd), sub("_[^_]+$", "", sp.n), las=2)
axis(side=2, at=seq(0, 20, 4))
tmp <- seq(0, 3e9, 1e9)
axis(side=4, at=tmp*tmp.gm, labels=tmp, line=0, las=2)
mtext("log 2 intron length", 2, line=3)

## for a range of thresholds:
small.r.l <- sapply(5:10, function(th){
    tapply( intron.s, gene.stats$sp, function(x){
        s <- unlist(x)
        sum( log2(s) < th ) / length(s) })
})

par(mfrow=c(2,3))
for(i in 1:ncol(small.r.l)){
    plot( log2(genome.sizes[sp.n]), (small.r.l[,i])[sp.n] )
    usr <- par("usr")
    text( log2(genome.sizes[sp.n]), (small.r.l[,i])[sp.n], sub("_[^_]+$", "", sp.n),
         pos=ifelse( log2(genome.sizes[sp.n]) - usr[1] > usr[2] - log2(genome.sizes[sp.n]),
                     2,4 ))
}

intron.mean <- tapply( intron.s, gene.stats$sp, function(x){
    mean( log2(unlist(x)) ) })

i <- 4
par(mfrow=c(1,2))
plot( log2(genome.sizes[sp.n]), (small.r.l[,i])[sp.n], xlab='log2 genome size', ylab='sum(s < 256) / n' )
usr <- par("usr")
text( log2(genome.sizes[sp.n]), (small.r.l[,i])[sp.n], sub("_[^_]+$", "", sp.n),
     pos=ifelse( log2(genome.sizes[sp.n]) - usr[1] > usr[2] - log2(genome.sizes[sp.n]),
                2,4 ), cex=0.7)
##
plot( log2(genome.sizes[sp.n]), intron.mean[sp.n], xlab='log2 genome size', ylab='mean intron size' )
usr <- par("usr")
text( log2(genome.sizes[sp.n]), intron.mean[sp.n], sub("_[^_]+$", "", sp.n),
     pos=ifelse( log2(genome.sizes[sp.n]) - usr[1] > usr[2] - log2(genome.sizes[sp.n]),
                2,4 ), cex=0.7)

l <- length(sp.n)
int.med.sp.mi <- vector(mode='list', length=(l^2 - l)/2)

k <- 0
for(i in 1:(length(sp.n)-1)){
    for(j in (i+1):length(sp.n)){
        k <- k + 1
        int.med.sp.mi[[k]] <- c('i'=i, 'j'=j,'sp1'=sp.n[i], 'sp2'=sp.n[j],
                                intron.med.mi.sp( int.s=intron.s, gene.st=gene.stats, sp1=sp.n[i], sp2=sp.n[j], numBins=20 ))
    }
}

### to have an example of plot
par(mfrow=c(1,1))
par(mar=c(5.1, 4.1, 4.1, 2.1))
for(i in 1:length(int.med.sp.mi)){
    plot.med.mi( int.med.sp.mi[[i]], tform=log, ps.count=1 )
    inpt <- readline("next: ")
}


## to visualise we can make a matrix and make some sort of reasonable plot
int.med.sp.mi.m <- matrix(0, nrow=length(sp.n), ncol=length(sp.n))
for(i in 1:length(int.med.sp.mi))
    int.med.sp.mi.m[ int.med.sp.mi[[i]]$i, int.med.sp.mi[[i]]$j ] = int.med.sp.mi[[i]]$mi


## but more control if I do:
y <- matrix( 1:length(sp.n), nrow=length(sp.n), ncol=length(sp.n), byrow=TRUE )
x <- matrix( 1:length(sp.n), nrow=length(sp.n), ncol=length(sp.n), byrow=FALSE )

cm <- 0.3  ## something to moderate the color
par(mar=c(1.1,1.1,1.1,1.1))
plot.new()
plot.window(xlim=c(-3,length(sp.n)), ylim=c(1, 5+length(sp.n)))
rect( x, y, x+1, y+1, col=hsvScale(int.med.sp.mi.m, val=(cm + int.med.sp.mi.m)/(cm + max(int.med.sp.mi.m))  ), border='grey', lwd=0.25 )
text( 0.95, y[1,] + 0.5, sub("_[^_]+$", "", sp.n), adj=c(1,0.5))
text( x[,1]+0.5, length(sp.n) + 1.05, sub("_[^_]+$", "", sp.n), adj=c(0,0.5), srt=90)


## convert to two different kinds of distance measures
int.med.sp.mi.d.cr <- matrix(0, nrow=length(sp.n), ncol=length(sp.n))
int.med.sp.mi.d.cl <- matrix(0, nrow=length(sp.n), ncol=length(sp.n))
for(i in 1:length(int.med.sp.mi)){
    cr <- 1 - (int.med.sp.mi[[i]]$mi / int.med.sp.mi[[i]]$mi.i)
    cl <- 1 - (int.med.sp.mi[[i]]$mi / max(int.med.sp.mi[[i]]$H1, int.med.sp.mi[[i]]$H2))
    r <- int.med.sp.mi[[i]]$i
    c <- int.med.sp.mi[[i]]$j
    int.med.sp.mi.d.cr[r,c] = cr
    int.med.sp.mi.d.cr[c,r] = cr
    int.med.sp.mi.d.cl[r,c] = cl
    int.med.sp.mi.d.cl[c,r] = cl
}

cm <- 0.3  ## something to moderate the color
par(mar=c(1.1,1.1,1.1,1.1))
plot.new()
plot.window(xlim=c(-3,length(sp.n)), ylim=c(1, 5+length(sp.n)))
rect( x, y, x+1, y+1, col=hsvScale(int.med.sp.mi.d.cr, val=0.8 * (int.med.sp.mi.d.cr)/(cm + max(int.med.sp.mi.d.cr))  ), border='grey', lwd=0.25 )
text( 0.95, y[1,] + 0.5, sub("_[^_]+$", "", sp.n), adj=c(1,0.5))
text( x[,1]+0.5, length(sp.n) + 1.05, sub("_[^_]+$", "", sp.n), adj=c(0,0.5), srt=90)


require("ape")

dimnames( int.med.sp.mi.d.cr ) <- list( sp.n, sp.n )
int.med.sp.mi.d.cr.nj <- nj(int.med.sp.mi.d.cr)

plot(int.med.sp.mi.d.cr.nj, "u")
plot(int.med.sp.mi.d.cr.nj, "fan")
plot(int.med.sp.mi.d.cr.nj, "radial")
plot(root(int.med.sp.mi.d.cr.nj, outgroup=sp.n[ grep("eptatretus", sp.n) ] )) ## this is probably the best way to do it..

plot(root(int.med.sp.mi.d.cr.nj, outgroup=sp.n[18]), 'radial') ## this is probably the best way to do it..

#### we can compare the intron distances to phylogenetic distances using the Ensembl
#### tree

sp.tree <- get.ensembl.taxonomy( db="ensembl_compara_98", root.name="Ensembl" )
rownames(sp.tree) <- sp.tree[,'node_id']


sp.names <- sub("_", " ", sp.n)
## side effects in R. Who would have through it..
substr(sp.names, 1, 1) <- toupper(substr(sp.names, 1, 1))

sum( sp.names %in% sp.tree$node_name ) ## 166
## we need to do a little modification ..
sp.names[ !(sp.names %in% sp.tree$node_name ) ]
## [1] "Canis familiaris" "Mus pahari"       "Mus caroli"       "Mus spretus"     
## [5] "Cebus capucinus"  "Gorilla gorilla"
## modify as follows
sp.names[ sp.names == 'Canis familiaris' ] <- 'Canis lupus'
sp.names[ sp.names == 'Mus pahari' ] <- 'Mus pahari strain PAHARI_EIJ'
sp.names[ sp.names == 'Mus caroli' ] <- 'Mus caroli strain CAROLI_EIJ'
sp.names[ sp.names == 'Mus spretus' ] <- 'Mus spretus strain SPRET/EiJ'
sp.names[ sp.names == 'Cebus capucinus' ] <- 'Cebus capucinus imitator'
sp.names[ sp.names == 'Gorilla gorilla' ] <- "Gorilla gorilla gorilla"

sum( sp.names %in% sp.tree$node_name ) ## 172

table(sapply( sp.names, function(x){ sum( sp.tree$node_name == x ) }))
## some are present more than once in the tree..
##   1   2  12  15 
## 169   1   1   1 

## create a vector in the same way that we calculated mutual information..
l <- length(sp.n)
ens.dist <- vector(mode='list', length=(l^2 - l)/2 )

k <- 0
for(i in 1:(length(sp.n)-1)){
    for(j in (i+1):length(sp.n)){
        k <- k + 1
        d <- get.ensembl.dist( sp.tree, sp.names[i], sp.names[j] )
        ens.dist[[k]] <- list('i'=i, 'j'=j,'sp1'=sp.n[i], 'sp2'=sp.n[j],
                           'd1'=d[1], 'd2'=d[2], 'd'=d[3] )
    }
}

## let us get the lineage for all of the nodes
sp.ncbi.tree <- get.ensembl.taxonomy( root.name='NCBI Taxonomy' )
sum( sp.names %in% sp.ncbi.tree$node_name ) ## 172
sp.lineages <- lapply( sp.names, function(x){ get.ensembl.lineage( sp.ncbi.tree, x ) } )
names(sp.lineages) <- sp.names 

saveRDS(sp.lineages, file='species_lineages.rds')

osteo.b <- as.logical(sapply( sp.lineages, function(x){ sum(x$node_name %in% "Osteoglossocephalai") } ))
mammal.b <- as.logical(sapply( sp.lineages, function(x){ sum(x$node_name %in% "Mammalia") } ))
aves.b <- as.logical(sapply( sp.lineages, function(x){ sum(x$node_name %in% "Aves") } ))

a <- sapply( int.med.sp.mi, function(x){ x$i })
b <- sapply( int.med.sp.mi, function(x){ x$j })

alpha <- 0.5
cols <- rep(rgb(0,0,0,alpha), length(x))
## this could be done better with bitwise operations, 1, 2, 4, 8
## for all, fish, mammal, bird
## but do it the stupid way for now
cols[ osteo.b[a] & osteo.b[b] ] <- rgb(1, 0, 0, alpha)
cols[ mammal.b[a] & mammal.b[b] ] <- rgb(0, 1, 0, alpha)
cols[ aves.b[a] & aves.b[b] ] <- rgb(0, 1, 0, alpha)
cols[ (osteo.b[a] & mammal.b[b]) | (osteo.b[b] & mammal.b[a])  ] <- rgb(1, 1, 0, alpha)
cols[ (osteo.b[a] & aves.b[b]) | (osteo.b[b] & aves.b[a])  ] <- rgb(1, 0, 1, alpha)
cols[ (aves.b[a] & mammal.b[b]) | (aves.b[b] & mammal.b[a])  ] <- rgb(0, 1, 1, alpha)

par('mar'=c(5.1, 4.1, 4.1, 2.1))
x <- sapply( int.med.sp.mi, function(x){ x$mi })
y <-  sapply( ens.dist, function(x){ x$d } )
plot( x, y, xlab='Mutual information', ylab='Ensembl distance',
     col=cols,
     panel.first=with(par(), rect(usr[1], usr[3], usr[2], usr[4], col=rgb(0.5, 0.5, 0.5))),
     lwd=2)

identify( x, y, sapply( ens.dist, function(x){ paste( x$sp1, x$sp2, sep='\n' ) }))

## the patterns are not as clear as before, but there is 


### Get N50 / N90 etc from the assembly golden paths
assemblies <- lapply( dbs, get.assembly )
saveRDS(assemblies, file='assemblies.rds')

## unfortunatly many of those are empty, so we will need to use
## the seq_region instead

## lots of 
seq.region <- lapply( dbs, get.chr.regions )

plot( log2(genome.sizes), log2( 1 + sapply(seq.region, function(x){ sum(x$length) })))
abline(0,1)
## there are a whole load of regions here that are not good

identify( log2(genome.sizes), log2(1 + sapply(seq.region, function(x){ sum(x$length) })),
         names(genome.sizes))
##
## oreochromis_niloticus has far too many seq_region and there is no

## we can get rid of most of the problem by:
b <- log2(genome.sizes) - log2( 1 + sapply(seq.region, function(x){ sum(x$length) })) < -1

seq.region[b] <- lapply( seq.region[b], function(x){ x[ !grepl("_", x$name), ] })

plot( log2(genome.sizes), log2( 1 + sapply(seq.region, function(x){ sum(x$length) })))
abline(0,1)

## oreochromis_niloticus is now left with far too many contigs.
with( seq.region[['oreochromis_niloticus']], {
    prefix <- unique( substr(name, 1, 5))
    y <- sapply( prefix, function(x){ sum( length[ grepl(x, name)] )})
    names(y) <- prefix
    y / genome.sizes['oreochromis_niloticus']
})

##        AERX0        GL831        GL832        AAQR0        GL874        GL875 
## 0.0086610705 0.9798405707 0.0114983588 0.0131277539 0.0097371163 0.0001921417 
##        GL873 
## 2.6939687255 

## so we can probably use GL831

b <- with(seq.region[['oreochromis_niloticus']], grep("GL831", name))
seq.region[['oreochromis_niloticus']] <- seq.region[['oreochromis_niloticus']][ b, ]


plot( log2(genome.sizes), log2( 1 + sapply(seq.region, function(x){ sum(x$length) })))
abline(0,1)

## this now looks more reasonable. We can try to export this.
saveRDS(seq.region, "seq_region.rds")

