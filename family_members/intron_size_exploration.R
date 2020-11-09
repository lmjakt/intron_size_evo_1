## This looks at intron sizes in a set of species from Ensembl.
## Here we will do:

## 1. Repeat some of the plots I created that I performed when looking at data
##    from recipocral blast for L. piscatorius. This was actually done to confirm
##    observations of bimodal intron size distributions made by Arseny. I also extended
##    the data set to include other species.
## 2. Extract family identifiers from the selected species and potentially make some selection
##    from these which will allow us to identify orthologous intron sequences. For example
##    such sequences are most easily identified when we have strict one-to-one orthology.
##    However, that in itself may be difficult given historical genome duplications
##    giving us not only paralogue creation, but also random loss at different times.
## 3. Output reasonable gene sets that should be aligned to each other. These identifiers
##    will then be used to extract sequences (probably using a perl script as I am not
##    sure of how this can be done efficiently in R) where we store the information
##    of exon intron boundaries, so that we can extract the positions of introns from
##    these alignments.
## 4. Again, probably in a perl script. Make use of the alignments in 3 to extract
##    orthologous intron sequnces and their lengths. Then we can first ask, whether
##    intron length is conserved. That is, are the set of introns shrinking in different
##    species of teleosts.

source("functions.R")
require(RMySQL)
## to clear all database connections
invisible( lapply( dbListConnections(dbDriver(drv="MySQL") ), dbDisconnect))

db.names <- c('ciona_intestinalis_core_97_3', 'danio_rerio_core_97_11', 'gasterosteus_aculeatus_core_97_1', 'lepisosteus_oculatus_core_97_1',
              'oryzias_latipes_core_97_1', 'takifugu_rubripes_core_97_5', 'tetraodon_nigroviridis_core_97_8',
              'poecilia_formosa_core_97_512', 'ursus_americanus_core_97_1', 'anolis_carolinensis_core_97_2',
              'lonchura_striata_domestica_core_97_1', 'callorhinchus_milii_core_97_613', 'electrophorus_electricus_core_97_2',
              'eptatretus_burgeri_core_97_32', 'neolamprologus_brichardi_core_97_1', 'amphilophus_citrinellus_core_97_5',
              'mus_musculus_core_97_38', 'homo_sapiens_core_97_38', 'saccharomyces_cerevisiae_core_97_4', 'oreochromis_niloticus_core_97_1',
              'mastacembelus_armatus_core_97_11', 'caenorhabditis_elegans_core_97_269', 'drosophila_melanogaster_core_97_7',
              'gallus_gallus_core_97_6', 'taeniopygia_guttata_core_97_1')

names(db.names) <- c('ciona', 'd.rerio', 'stickleback', 'gar', 'medaka', 'fugu', 't.odon',
                     'a.molly', 'bl.bear', 'a.lizard', 'b.finch', 'e.shark', 'e.eel', 'hagfish', 'l.cichlid', 'm.cichlid',
                     'mus', 'human', 'yeast', 'tilapia', 'zz.eel', 'c.elegans', 'f.fly', 'chicken', 'z.finch')
sp.names <- sub("_core_97_\\d+$", "", db.names)

## we can't connect to all of these at the same time, as it exceeds the number of database connections that
## I am allowed. This means that we have a lot of connecting and disconnecting time, but it is not that
## important.

## for figure proportions
a5.w <- 5.845
a5.h <- 8.27
im.m <- 1.5

## lets get genome stats:
genome.stats <- lapply( db.names, get.genome.stats )
genome.sizes <- sapply( genome.stats, function(x){ x[ x$statistic == 'ref_length', 'value']  })
size.o <- order(genome.sizes)  ## use this later..

## then let us get the exon information. This will take a bit longer.
exon.info <- lapply( db.names, get.ens.exon.info )

## then let us get the lengths from these..
intron.lengths <- lapply( exon.info, extract.intron.lengths )

## then let us visualise the data in a reasonable manner
intron.lengths.all.h <- hist( log2( unlist( sapply( intron.lengths, function(x){ x[,'size'] }))), breaks=100 )

intron.lengths.h <- lapply( intron.lengths, function(x){ hist( log2(x[,'size']), breaks=intron.lengths.all.h$breaks )})
intron.lengths.h.1 <- lapply( intron.lengths, function(x){ hist( log2(x[x[,'rank'] == 1,'size']), breaks=intron.lengths.all.h$breaks )})
intron.lengths.h.2 <- lapply( intron.lengths, function(x){ hist( log2(x[x[,'rank'] > 1, 'size']), breaks=intron.lengths.all.h$breaks )})
intron.lengths.h.3 <- lapply( intron.lengths, function(x){ hist( log2(x[x[,'irank'] == 1, 'size']), breaks=intron.lengths.all.h$breaks )})

pdf('first_intron_size_hists.pdf', width=a5.w * im.m, height=a5.h * im.m, title='Intron size distributions')
plot.intron.hist( intron.lengths.h[size.o], intron.lengths.h.1[size.o], intron.lengths.h.2[size.o], lab.cex=1,
                 v.lines=intron.lengths.h[[1]]$mids[ seq(10, 100, 10) ], v.lin.lwd=0.8 )
dev.off()

## apart from the cichlids we do not see much of a difference for the terminal exon. So we can leave this out
## here.
plot.intron.hist( intron.lengths.h[size.o], intron.lengths.h.1[size.o], intron.lengths.h.3[size.o], lab.cex=1,
                 v.lines=intron.lengths.h[[1]]$mids[ seq(10, 100, 10) ], v.lin.lwd=0.8 )

## Get genes that map to non-variant scaffolds or chromosomes. This in order
## to avoid overcounting family members
ref.genes <- lapply( db.names, get.all.ref.genes, max.coord.rank=2 )


#############################################################################################################
### Get gene family members for these set of species from Compara ###########################################

compara <- ensembl.connect("ensembl_compara_97")

family.members <- lapply(sp.names, get.family.members, db=compara)
families <- dbGetQuery( compara, "select * from family" )

## let us remove alternate genes
sapply(1:length(family.members), function(i){
    c(length(unique(family.members[[i]]$gene.id)), sum( unique(family.members[[i]]$gene.id) %in% ref.genes[[i]]$stable_id) ) })

##       [,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10] [,11] [,12]
## [1,] 10307 28450 19250 17622 21811 19306 18925 22330 19414 17718 14813 16963
## [2,] 10307 23868 17045 17622 21811 19306 18925 22330 19414 16762 14813 16963
##      [,13] [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24]
## [1,] 20922 13465 21528 21661 21972 23620  5363 20779 21729 15276 13117 16225
## [2,] 20922 13465 21528 21661 21647 19973  5363 20779 21729 15276 13117 16225
##      [,25]
## [1,] 16365
## [2,] 16365

for(i in 1:length(family.members)){
    family.members[[i]] <- family.members[[i]][ family.members[[i]][,'gene.id'] %in% ref.genes[[i]][,'stable_id'], ]
}

## let us get a table of counts for each species..
family.counts <- sapply( family.members, function(x){
    table(c(families$stable_id, x$family.id[!duplicated( x$gene.id )])) - 1 })

dim(family.counts)
## [1] 32296    25
family.counts.h <- hist( rowSums(family.counts == 1) )
family.counts.q <- sort( rowSums(family.counts == 1) )

## that is a bit hopeless. I suppose that we should only look at vertebrate species
vert.b <- genome.sizes > 300e6
family.counts.v.h <- hist( rowSums(family.counts[,vert.b] == 1), breaks=0:22-0.5 ) ## only 2 families with one gene in all
family.counts.v.h.2 <- hist( rowSums(family.counts[,vert.b] > 0 & family.counts[,vert.b] < 3), breaks=0:22-0.5 ) ##
family.counts.v.h.4 <- hist( rowSums(family.counts[,vert.b] > 0 & family.counts[,vert.b] < 5), breaks=0:22-0.5 )

## 
table( rowSums(family.counts[,vert.b] > 0 & family.counts[,vert.b] < 5) )
##     0     1     2     3     4     5     6     7     8     9    10    11    12 
## 13347  1421  1041  1466   553   400   365   453   375   376   386   490   544 
##    13    14    15    16    17    18    19    20    21 
##   347   328   345   463   639   915  1566  2832  3644 

table( rowSums(family.counts[,vert.b] > 0 & family.counts[,vert.b] == 1) )
##     0     1     2     3     4     5     6     7     8     9    10    11    12 
## 14245  1816  1271  1479   698   564   605   688   722   796   796   717   564 
##    13    14    15    16    17    18    19    20    21 
##   415   437   491   621   863  1168  1416  1234   690 

## if we take orthologues from at least 15 species we get:
## this 
491 + 621 + 863 + 1168 + 1416 + 1234 + 690
## [1] 6483

## 
table( rowSums(family.counts[,vert.b] > 0 & family.counts[,vert.b] < 5) )
##     0     1     2     3     4     5     6     7     8     9    10    11    12 
## 13347  1421  1041  1466   553   400   365   453   375   376   386   490   544 
##    13    14    15    16    17    18    19    20    21 
##   347   328   345   463   639   915  1566  2832  3644 

## We can start the analysis with strict single member families, and at least
## 15 species for each one. This gets us a total 6483 genes, and presumably
## rather more 

## this is incredibly inefficient, but never mind.
fam.1.15.members <- lapply( rownames(family.counts)[ rowSums(family.counts[,vert.b] > 0 & family.counts[,vert.b] == 1) >= 15 ], function(fam){
    lapply( family.members[ vert.b ], function(x){
        unique( x$gene.id[ x$family.id == fam ])
    })
})

names(fam.1.15.members) <- rownames(family.counts)[ rowSums(family.counts[,vert.b] > 0 & family.counts[,vert.b] == 1) >= 15 ]

## that actually allows for some gene duplications within a species; but most will be single genes. Let us output this list
## as a table
fam.1.15.tbl <- t(sapply( fam.1.15.members, function(x){
    sapply(x, paste, collapse=',') }))

colnames(fam.1.15.tbl) <- db.names[ vert.b ]

write.table( fam.1.15.tbl, "vertebrate_family_members_1_15.txt", quote=FALSE, sep="\t" )
