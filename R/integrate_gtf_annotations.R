## Andrea Hita
## OCT-2020

## Import libraries (these packages are required for the gtf integration)
# library(Hmisc)
# library(parallel)
# library(plyr)
# library(dplyr)
# library(rtracklayer)

## Paths definition
##db.dir <- "~"
    
## Script parameters
##nc <- 20

## Define analysis flags
annotHuman <- FALSE
annotMouse <- FALSE
annotArabidopsis <- FALSE
annotNematode <- FALSE

## BIOCATEGORIES DEFINITION
## ----------------------------------------------------------------------------

## Define biotype categories and merge it with .gtf
define_bcats <- function(){
    
    bcats <- data.frame(rbind(

    ## Immunoglobulin (Ig) genes
    cbind("Protein_coding", "IG", c("IG_C_gene","IG_D_gene","IG_J_gene",
                                    "IG_V_gene","IG_LV_gene")),

    ## Innactivated Immunoglobulin pseudogenes
    cbind("Protein_coding", "IG_p", c("IG_pseudogene",
                                      "IG_C_pseudogene","IG_D_pseudogene",
                                      "IG_J_pseudogene","IG_V_pseudogene")),
        
    ## T-cell receptor genes
    cbind("Protein_coding", "TR", c("TR_C_gene","TR_D_gene","TR_J_gene","TR_V_gene")),

    ## Innactivated T-cell receptor pseudogenes
    cbind("Protein_coding", "TR_p", c("TR_J_pseudogene", "TR_V_pseudogene")),

    ## Protein coding genes
    cbind("Protein_coding", "Prot","protein_coding"),

    ## Small non-coding rna
    cbind("Short_non_coding", "sNC", c("snRNA", "snoRNA", "scaRNA", "pre_miRNA", "miRNA", "Y_RNA", "siRNA",
                                       "vault_RNA", "ribozyme", "misc_RNA", "piRNA", "scRNA", "sRNA")),

    ## Long non-coding
    cbind("Long_non_coding", "lNC", c("lncRNA","lincRNA","macro_lncRNA","ncRNA",
                                      "antisense","sense_intronic","sense_overlapping",
                                      "processed_transcript","3prime_overlapping_ncRNA",
                                      "bidirectional_promoter_lncRNA","nontranslating_CDS",
                                      "antisense_RNA")),    

    ## Ribosomal rna genes
    cbind("rRNA", "rRNA", "rRNA"),

    ## Ribosomal rna pseudogenes
    cbind("rRNA", "rRNA_p", c("rRNA_pseudogene")),

    ## Mithocondrial rna
    cbind("rRNA", "Mt_rRNA", "Mt_rRNA"),

    ## Ribosomal rna genes
    cbind("tRNA", "tRNA", c("tRNA","tRF","tRF5","tRF3")),

    ## Mithocondrial rna
    cbind("tRNA", "Mt_tRNA", "Mt_tRNA"),

    ## Other pseudogenes
    cbind("Long_pseudogenes", "Other_p",
          c("pseudogene", "polymorphic_pseudogene",
            "unitary_pseudogene",
            "unprocessed_pseudogene","processed_pseudogene",
            "transcribed_unprocessed_pseudogene",
            "transcribed_processed_pseudogene",
            "transcribed_unitary_pseudogene",
            "translated_unprocessed_pseudogene",
            "translated_processed_pseudogene")),

    cbind("Hybrid", "TEC", "TEC"),
    cbind("Hybrid", "Hybrid", "Hybrid")))

    names(bcats) <- c("biogroup", "biocat", "biotype")
    return(bcats)}

bcats <- define_bcats()
##write.csv(bcats, file = paste0(db.dir,"bcats.csv"), row.names = FALSE)

## HUMAN ANNOTATIONS EXTRACTION
## ----------------------------------------------------------------------------

if(annotHuman){

    ## ---- Download annotations
    download.file(url = "ftp://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.103.gtf.gz",
                  destfile = paste0(db.dir,"ensembl_Homo_sapiens.GRCh38.103.gtf.gz"))
    download.file(url = "ftp://mirbase.org/pub/mirbase/CURRENT/genomes/hsa.gff3",
                  destfile = paste0(db.dir,"mirbase_hsa.gff3"))

    ## ---- Import raw annotations
    print("Importing annotations...")
    gtf.ens <- rtracklayer::import(paste0(db.dir,"ensembl_Homo_sapiens.GRCh38.103.gtf.gz"))
    gff.dashr <- rtracklayer::import(paste0(db.dir,"dashr.v2.sncRNA.annotation.hg38.gff"))
    gff.mir <- rtracklayer::import(paste0(db.dir,"mirbase_hsa.gff3"))
    
    ## ---- 1. Curate Ensembl annotations

    ## Remove annotations outside chr and MT sequences
    gtf <- gtf.ens[as.character(gtf.ens@seqnames) %in% levels(gtf.ens@seqnames)[1:25],]
    gtf@seqnames@values <- droplevels(gtf@seqnames@values)
    
    ## Convert less-represented categories to miscRNA
    gtf@elementMetadata$gene_biotype[gtf@elementMetadata$gene_biotype == "scRNA"] <- "misc_RNA"
    gtf@elementMetadata$gene_biotype[gtf@elementMetadata$gene_biotype == "sRNA"] <- "misc_RNA"
    gtf@elementMetadata$gene_biotype[gtf@elementMetadata$gene_biotype == "vault_RNA"] <- "misc_RNA"

    gtf@elementMetadata$transcript_biotype[gtf@elementMetadata$transcript_biotype == "scRNA"] <- "misc_RNA"
    gtf@elementMetadata$transcript_biotype[gtf@elementMetadata$transcript_biotype == "sRNA"] <- "misc_RNA"
    gtf@elementMetadata$transcript_biotype[gtf@elementMetadata$transcript_biotype == "vault_RNA"] <- "misc_RNA"
    
    ## Exclude miRNA annotations
    gtf <- gtf[gtf@elementMetadata$gene_biotype != "miRNA",]

    ## ...exclude 2 miRNA annotations wrongly annotated as lncRNA
    gtf <- gtf[gtf@elementMetadata$gene_id != "ENSG00000272920",]
    gtf <- gtf[gtf@elementMetadata$gene_id != "ENSG00000266919",]
    
    ## ...exclude 2 lncRNA annotations that get rRNA
    gtf <- gtf[gtf@elementMetadata$gene_id != "ENSG00000278996",]
    gtf <- gtf[gtf@elementMetadata$gene_id != "ENSG00000280441",]
    
    ## ---- 2. Re-format pre-miRNA miRbase gff annotations to .gtf

    print("Integrating miRNA from miRbase...")
    
    ## Relabel sequence names
    gff.mir@seqnames@values <- mapvalues(gff.mir@seqnames@values, from = levels(gff.mir@seqnames),
                                         to = c(as.character(c(1:22)),"X","Y"))

    gff.mir.p <- gff.mir[gff.mir@elementMetadata$type == "miRNA_primary_transcript"]
    gff.mir.m <- gff.mir[gff.mir@elementMetadata$type == "miRNA"]

    ## --------------------

    ## pre-miRNA annotations
    
    ## Modify gene_name to avoid duplicates
    gid.dup <-  gff.mir.p@elementMetadata$Name[which(duplicated(gff.mir.p@elementMetadata$Name) == TRUE)]
    gidx.d <- which(gff.mir.p@elementMetadata$Name %in% gid.dup)
    gff.mir.p[gidx.d,]

    tmp <- gff.mir.p
    gff.mir.p@elementMetadata$Name[gidx.d[1:4]] <- paste(gff.mir.p@elementMetadata$Name[gidx.d[1:4]],c(1:4),sep="-")
    gff.mir.p <- gff.mir.p[-gidx.d[c(6,8)],]

    gidx.d.m <- which(gff.mir.m@elementMetadata$Derives_from %in% tmp@elementMetadata$ID[gidx.d])
    gff.mir.m[gidx.d.m[7:12],]
    gff.mir.m@elementMetadata$Derives_from[gidx.d.m[c(1:8,12)]] <- paste0(
        gff.mir.m@elementMetadata$Derives_from[gidx.d.m[c(1:8,11)]],c("","","_2","_2","_3","_3","_4","_4","_2"))
    gff.mir.m <- gff.mir.m[-gidx.d.m[10:12],]
    
    ## Create an Ensemble like gtf empty templage
    gtf.template <- gtf.ens[1,]
    gtf.template@elementMetadata[,-2] <- NA
    gtf.template@elementMetadata$gene_biotype <- "pre_miRNA"
    gtf.template@elementMetadata$transcript_biotype <- NA
    gtf.template@elementMetadata$source <- "mirbase"
    
    ## pre-miRNA .gtf generation loop
    out <- mclapply(c(1:length(gff.mir.p@seqnames)), function(i) {

        ## Initialize .gtf template
        gtf.mir.p <- gtf.template

        ## Define annot coordinates
        gtf.mir.p@ranges <- gff.mir.p@ranges[i]
        gtf.mir.p@seqnames <- gff.mir.p@seqnames[i]
        gtf.mir.p@strand <- gff.mir.p@strand[i]

        ## Define metadata 
        gtf.mir.p@elementMetadata$gene_id <- gff.mir.p@elementMetadata$ID[i]
        gtf.mir.p@elementMetadata$gene_name <- gff.mir.p@elementMetadata$Name[i]

        return(gtf.mir.p)}, mc.cores = nc)

    ## Concatenate gtf.mir.p 
    gtf.mir.p <- gtf.template
    slot(gtf.mir.p,"seqnames") <- do.call("c", mclapply(out, function(x) x@seqnames, mc.cores = nc))
    slot(gtf.mir.p,"strand") <- do.call("c", mclapply(out, function(x) x@strand, mc.cores = nc))
    slot(gtf.mir.p,"ranges") <- do.call("c", mclapply(out, function(x) x@ranges, mc.cores = nc))
    slot(gtf.mir.p,"elementMetadata")  <- do.call("rbind", mclapply(out, function(x)
        x@elementMetadata, mc.cores = nc))

    ## --------------------
    
    ## Mature miRNA annotations
        
    ## Create an Ensemble like gtf empty templage
    gtf.template <- gtf.ens[1:2,]
    gtf.template@elementMetadata[,] <- NA
    gtf.template@elementMetadata$gene_biotype <- "pre_miRNA"
    gtf.template@elementMetadata$transcript_biotype <- "miRNA"
    gtf.template@elementMetadata$source <- "mirbase"
    gtf.template@elementMetadata$type <- c("transcript","exon")

    ## mature-miRNA .gtf generation loop
    out <- mclapply(c(1:length(gff.mir.m@seqnames)), function(i) {

        ## Initialize .gtf template
        gtf.mir.m <- gtf.template

        ## Define annot coordinates
        gtf.mir.m@ranges[c(1,2)] <- gff.mir.m@ranges[i]
        gtf.mir.m@seqnames <- c(gff.mir.m@seqnames[i], gff.mir.m@seqnames[i])
        gtf.mir.m@strand <- c(gff.mir.m@strand[i], gff.mir.m@strand[i])

        ## Define metadata
        gtf.mir.m@elementMetadata$gene_id <- gff.mir.m@elementMetadata$Derives_from[i]
        gtf.mir.m@elementMetadata$transcript_id <- gff.mir.m@elementMetadata$ID[i]
        gtf.mir.m@elementMetadata$transcript_name <- gff.mir.m@elementMetadata$Name[i]
        gtf.mir.m@elementMetadata$exon_id[2] <- paste0(gff.mir.m@elementMetadata$ID[i],".e")
        gtf.mir.m@elementMetadata$gene_name <-
            gff.mir.p@elementMetadata$Name[gff.mir.p@elementMetadata$ID ==
                                           gff.mir.m@elementMetadata$Derives_from[i]]

        return(gtf.mir.m)}, mc.cores = nc)

    ## Concatenate gtf.mir.m 
    gtf.mir.m <- gtf.template
    slot(gtf.mir.m,"seqnames") <- do.call("c", mclapply(out, function(x) x@seqnames, mc.cores = nc))
    slot(gtf.mir.m,"strand") <- do.call("c", mclapply(out, function(x) x@strand, mc.cores = nc))
    slot(gtf.mir.m,"ranges") <- do.call("c", mclapply(out, function(x) x@ranges, mc.cores = nc))
    slot(gtf.mir.m,"elementMetadata")  <- do.call("rbind", mclapply(out, function(x)
        x@elementMetadata, mc.cores = nc))

    ## Modify transcript_name to avoid duplicates
    gtf.mir.m.t <- subset(gtf.mir.m, type == 'transcript')
    tid.dup <- gtf.mir.m.t@elementMetadata$transcript_name[
                 which(duplicated(gtf.mir.m.t@elementMetadata$transcript_name) == TRUE)]
    tidx.d <- which(gtf.mir.m@elementMetadata$transcript_name %in% tid.dup)

    suff <- sapply(tidx.d, function(i)
        sub(as.character(sub("-3p","",sub("-5p","", sub(
              "R","r",as.character(gtf.mir.m@elementMetadata$transcript_name[i]))))),
            "", gtf.mir.m@elementMetadata$gene_name[i]))
    suff[599:600] <- ""
    gtf.mir.m@elementMetadata$transcript_name[tidx.d] <- paste0(
        gtf.mir.m@elementMetadata$transcript_name[tidx.d],suff)

    ## Manual correction of miss-annotation
    gtf.mir.m@elementMetadata$transcript_name[gtf.mir.m@elementMetadata$transcript_id ==
                                              "MIMAT0022737"] <- "hsa-miR-550b-1-5p"
    
    ## --------------------
    ## ---- 3. Re-format piRNA DASHR .gff annotations to .gtf

    print("Integrating piRNA from DASHR...")
    
    ##  Get piRNA from DASHR database 
    gff.pi <- gff.dashr[gff.dashr@elementMetadata$type == "piRNA",]
    gff.pi@seqnames@values <- droplevels(gff.pi@seqnames@values)
    gff.pi@seqnames@values <- mapvalues(gff.pi@seqnames@values, from = levels(gff.pi@seqnames),
                                        to = c(as.character(c(1:22)),"X","Y","MT"))
    for(i in unique(gff.pi@elementMetadata$ID[duplicated(gff.pi@elementMetadata$ID)])){
       idx <- which(gff.pi@elementMetadata$ID == i)
       gff.pi@elementMetadata$ID[idx] <- paste0(gff.pi@elementMetadata$ID[idx],"-",
                                                c(1:length(idx)))}

    ## Create an Ensemble like gtf empty templage
    gtf.template <- gtf.ens[1:3,]
    gtf.template@elementMetadata[1:3,-2] <- NA
    gtf.template@elementMetadata$gene_biotype <- rep("piRNA", 3)
    gtf.template@elementMetadata$transcript_biotype <- c(NA, rep("piRNA", 2))
    gtf.template@elementMetadata$source <- "DASHR"

    ## piRNA .gtf generation loop
    out <- mclapply(c(1:length(gff.pi@seqnames)), function(i) {

        ## Initialize .gtf template
        gtf.pi <- gtf.template

        ## Define annot coordinates
        gtf.pi@ranges[c(1:3)] <- gff.pi@ranges[i]
        gtf.pi@seqnames <- c(gff.pi@seqnames[i], gff.pi@seqnames[i], gff.pi@seqnames[i])
        gtf.pi@strand <- c(gff.pi@strand[i], gff.pi@strand[i], gff.pi@strand[i])

        ## Define metadata 
        ##gtf.pi@elementMetadata$gene_id[1:3] <- gff.pi@elementMetadata$ncbiAccession[i]
        ##gtf.pi@elementMetadata$transcript_id[2:3] <- paste0(gff.pi@elementMetadata$ncbiAccession[i],".t")
        ##gtf.pi@elementMetadata$exon_id[3] <- paste0(gff.pi@elementMetadata$ncbiAccession[i],".e")
        gtf.pi@elementMetadata$gene_id[1:3] <- sub("piR-","PIRNA",gff.pi@elementMetadata$ID[i])
        gtf.pi@elementMetadata$transcript_id[2:3] <- paste0(sub("piR-","PIRNA",gff.pi@elementMetadata$ID[i]),".t")
        gtf.pi@elementMetadata$exon_id[3] <- paste0(sub("piR-","PIRNA",gff.pi@elementMetadata$ID[i]),".e")
        gtf.pi@elementMetadata$gene_name[1:3] <- gtf.pi@elementMetadata$transcript_name[2:3] <-
            sub("piR-","PIRNA",gff.pi@elementMetadata$ID[i])

        return(gtf.pi)}, mc.cores = nc)

    ## Concatenate gtf.pi
    gtf.pi <- gtf.template
    slot(gtf.pi,"seqnames") <- do.call("c", mclapply(out, function(x) x@seqnames, mc.cores = nc))
    slot(gtf.pi,"strand") <- do.call("c", mclapply(out, function(x) x@strand, mc.cores = nc))
    slot(gtf.pi,"ranges") <- do.call("c", mclapply(out, function(x) x@ranges, mc.cores = nc))
    slot(gtf.pi,"elementMetadata")  <- do.call(
        "rbind", mclapply(out, function(x) x@elementMetadata, mc.cores = nc))

    ## --------------------
    ## ---- 4. Re-format tRF DASHR .gff annotations to .gtf

    print("Integrating tRF from DASHR...")
    
    ##  Get tRF from DASHR database 
    gff.trf <- gff.dashr[gff.dashr@elementMetadata$type %in% c("tRF5","tRNA","tRF3") &
                         gff.dashr@seqnames %in%
                         paste0('chr',c(as.character(c(1:22)),"X","Y","MT")),]
    gff.trf@seqnames@values <- droplevels(gff.trf@seqnames@values)
    gff.trf@seqnames@values <- mapvalues(gff.trf@seqnames@values, from = levels(gff.trf@seqnames),
                                         to = c(as.character(c(1:22)),"X"))
    
    ## Create an Ensemble like gtf empty templage
    gtf.template <- gtf.ens[1:7,]
    gtf.template@elementMetadata[,] <- NA
    gtf.template@elementMetadata$gene_biotype <- "tRNA"
    gtf.template@elementMetadata$transcript_biotype <- NA
    gtf.template@elementMetadata$source <- "DASHR"
    gtf.template@elementMetadata$type <- c("gene","transcript","exon","transcript","exon","transcript","exon")

    trfgenes <-  sub("-tRF3", "", sub("-tRF5", "", gff.trf@elementMetadata$ID))
    
    ## tRF .gtf generation loop
    out <- mclapply(unique(trfgenes), function(gi) {

        tt <- gff.trf[trfgenes %in% gi,]

        ## Initialize .gtf template
        gtf.trf <- gtf.template

        ## Define annot coordinates
        gtf.trf@ranges[1] <- IRanges(start = min(tt@ranges@start), end = max(tt@ranges@start + tt@ranges@width)-1)
        gtf.trf@ranges[2:3] <- tt@ranges[tt@elementMetadata$type == "tRF5"]
        gtf.trf@ranges[4:5] <- tt@ranges[tt@elementMetadata$type == "tRNA"]
        gtf.trf@ranges[6:7] <- tt@ranges[tt@elementMetadata$type == "tRF3"]
        gtf.trf@seqnames[1:7] <- c(tt@seqnames[1],tt@seqnames[1],tt@seqnames[1],tt@seqnames[1],
                                     tt@seqnames[1],tt@seqnames[1],tt@seqnames[1])
        gtf.trf@strand <- c(tt@strand[1],tt@strand[1],tt@strand[1],tt@strand[1],tt@strand[1],tt@strand[1],tt@strand[1])

         ## Define metadata
        gtf.trf@elementMetadata$transcript_biotype[2:3] <- "tRNA"
        gtf.trf@elementMetadata$transcript_biotype[4:5] <- "tRNA"
        gtf.trf@elementMetadata$transcript_biotype[6:7] <- "tRNA"
        gtf.trf@elementMetadata$gene_id[1:7] <- gtf.trf@elementMetadata$gene_name[1:7] <- gi
        gtf.trf@elementMetadata$transcript_id[2:3] <- gtf.trf@elementMetadata$transcript_name[2:3] <-
             tt@elementMetadata$ID[tt@elementMetadata$type == "tRF5"]
        gtf.trf@elementMetadata$transcript_id[4:5] <- gtf.trf@elementMetadata$transcript_name[4:5] <-
            tt@elementMetadata$ID[tt@elementMetadata$type == "tRNA"]
        gtf.trf@elementMetadata$transcript_id[6:7] <- gtf.trf@elementMetadata$transcript_name[6:7] <-
            tt@elementMetadata$ID[tt@elementMetadata$type == "tRF3"]
        gtf.trf@elementMetadata$exon_id[3] <-  paste0(gff.trf@elementMetadata$transcript_id[3],".e")
        gtf.trf@elementMetadata$exon_id[5] <-  paste0(gff.trf@elementMetadata$transcript_id[5],".e")
        gtf.trf@elementMetadata$exon_id[7] <-  paste0(gff.trf@elementMetadata$transcript_id[7],".e")

        return(gtf.trf)}, mc.cores = nc)

    ## Concatenate gtf.trf
    gtf.trf <- gtf.template
    slot(gtf.trf,"seqnames") <- do.call("c", mclapply(out, function(x) x@seqnames, mc.cores = nc))
    slot(gtf.trf,"strand") <- do.call("c", mclapply(out, function(x) x@strand, mc.cores = nc))
    slot(gtf.trf,"ranges") <- do.call("c", mclapply(out, function(x) x@ranges, mc.cores = nc))
    slot(gtf.trf,"elementMetadata")  <- do.call(
        "rbind", mclapply(out, function(x) x@elementMetadata, mc.cores = nc))
    
    ## --------------------
    ## ---- 5. Merge all .gtf annotations

    print("Exporting integrated annotations...")
    
    gtf.tmp <- gtf
    gtf.tmp@elementMetadata <- rbind(gtf.tmp@elementMetadata, gtf.pi@elementMetadata, gtf.trf@elementMetadata,
                                     gtf.mir.p@elementMetadata, gtf.mir.m@elementMetadata)
    gtf.tmp@seqnames <- c(gtf.tmp@seqnames, gtf.pi@seqnames, gtf.trf@seqnames, gtf.mir.p@seqnames, gtf.mir.m@seqnames)
    gtf.tmp@strand <- c(gtf.tmp@strand, gtf.pi@strand, gtf.trf@strand, gtf.mir.p@strand, gtf.mir.m@strand)
    gtf.tmp@ranges <- c(gtf.tmp@ranges, gtf.pi@ranges, gtf.trf@ranges, gtf.mir.p@ranges, gtf.mir.m@ranges)
    gtf <- gtf.tmp

    ## Export annotations
    gtf@elementMetadata$transcript_name <- sub('-201','',gtf@elementMetadata$transcript_name)
    rtracklayer::export(gtf, paste0(db.dir,'/Homo_sapiens.GRCh38.custom.gtf'))
    print("Human custom gtf successfully generated!")
}


## MOUSE ANNOTATIONS EXTRACTION
## ----------------------------------------------------------------------------

if(annotMouse){

    ## ---- Download annotations
    download.file(url = "ftp://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz",
                  destfile = paste0(db.dir,"ensembl_Mus_musculus.GRCm38.102.gtf.gz"))
    download.file(url = "ftp://mirbase.org/pub/mirbase/CURRENT/genomes/mmu.gff3",
                  destfile = paste0(db.dir,"mirbase_mmu.gff3"))
    download.file(url = "ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/16.0/genome_coordinates/gff3/mus_musculus.GRCm38.gff3.gz",
                  destfile = paste0(db.dir,"rnacentral_mus_musculus.GRCm38.gff3.gz"))
                  
    ## ---- Import raw annotations
    print("Importing annotations...")
    gtf.ens <- rtracklayer::import(paste0(db.dir,"ensembl_Mus_musculus.GRCm38.102.gtf.gz"))
    gff.mir <- rtracklayer::import(paste0(db.dir,"mirbase_mmu.gff3"))
    gff.rcent <- rtracklayer::import(paste0(db.dir,"rnacentral_mus_musculus.GRCm38.gff3.gz"))
        
     ## ---- 1. Curate Ensembl annotations
    
    ## Change underrpresented srna to miscRNA (only 1 and 2 of each respectively)
    gtf <-  gtf.ens[as.character(gtf.ens@seqnames) %in% levels(gtf.ens@seqnames)[1:22],]
    gtf@seqnames@values <- droplevels(gtf@seqnames@values)
    gtf@elementMetadata$gene_biotype[gtf@elementMetadata$gene_biotype == "scRNA"] <- "misc_RNA"
    gtf@elementMetadata$gene_biotype[gtf@elementMetadata$gene_biotype == "sRNA"] <- "misc_RNA"

    ## Exclude miRNA annotations
    gtf <- gtf[gtf@elementMetadata$gene_biotype != "miRNA",]

    ## ---- 2. Re-format pre-miRNA miRbase gff annotations to .gtf

    print("Integrating miRNA from miRbase...")
    
    ## Relabel sequence names
    gff.mir@seqnames@values <- mapvalues(gff.mir@seqnames@values, from = levels(gff.mir@seqnames),
                                         to = c(as.character(c(1:19)),"X"))
    levels(gff.mir@seqnames)[length(levels(gff.mir@seqnames))+1] <- "Y"

    gff.mir.p <- gff.mir[gff.mir@elementMetadata$type == "miRNA_primary_transcript"]
    gff.mir.m <- gff.mir[gff.mir@elementMetadata$type == "miRNA"]

    ## --------------------

    ## pre-miRNA annotations

    ## Create an Ensemble like gtf empty templage
    gtf.template <- gtf.ens[1,]
    gtf.template@elementMetadata[,-2] <- NA
    gtf.template@elementMetadata$gene_biotype <- "pre_miRNA"
    gtf.template@elementMetadata$transcript_biotype <- NA
    gtf.template@elementMetadata$source <- "mirbase"
    
    ## pre-miRNA .gtf generation loop
    out <- mclapply(c(1:length(gff.mir.p@seqnames)), function(i) {

        ## Initialize .gtf template
        gtf.mir.p <- gtf.template

        ## Define annot coordinates
        gtf.mir.p@ranges <- gff.mir.p@ranges[i]
        gtf.mir.p@seqnames <- gff.mir.p@seqnames[i]
        gtf.mir.p@strand <- gff.mir.p@strand[i]

        ## Define metadata 
        gtf.mir.p@elementMetadata$gene_id <- gff.mir.p@elementMetadata$ID[i]
        gtf.mir.p@elementMetadata$gene_name <- gff.mir.p@elementMetadata$Name[i]

        return(gtf.mir.p)}, mc.cores = nc)

    ## Concatenate gtf.mir.p 
    gtf.mir.p <- gtf.template
    slot(gtf.mir.p,"seqnames") <- do.call("c", mclapply(out, function(x) x@seqnames, mc.cores = nc))
    slot(gtf.mir.p,"strand") <- do.call("c", mclapply(out, function(x) x@strand, mc.cores = nc))
    slot(gtf.mir.p,"ranges") <- do.call("c", mclapply(out, function(x) x@ranges, mc.cores = nc))
    slot(gtf.mir.p,"elementMetadata")  <- do.call("rbind", mclapply(out, function(x) x@elementMetadata, mc.cores = nc))

    ## --------------------

    ## Mature miRNA annotations
        
    ## Create an Ensemble like gtf empty templage
    gtf.template <- gtf.ens[1:2,]
    gtf.template@elementMetadata[,] <- NA
    gtf.template@elementMetadata$gene_biotype <- "pre_miRNA"
    gtf.template@elementMetadata$transcript_biotype <- "miRNA"
    gtf.template@elementMetadata$source <- "mirbase"
    gtf.template@elementMetadata$type <- c("transcript","exon")

    ## mature-miRNA .gtf generation loop
    out <- mclapply(c(1:length(gff.mir.m@seqnames)), function(i) {

        ## Initialize .gtf template
        gtf.mir.m <- gtf.template

        ## Define annot coordinates
        gtf.mir.m@ranges[c(1,2)] <- gff.mir.m@ranges[i]
        gtf.mir.m@seqnames <- c(gff.mir.m@seqnames[i], gff.mir.m@seqnames[i])
        gtf.mir.m@strand <- c(gff.mir.m@strand[i], gff.mir.m@strand[i])

        ## Define metadata
        gtf.mir.m@elementMetadata$gene_id <- gff.mir.m@elementMetadata$Derives_from[i]
        gtf.mir.m@elementMetadata$transcript_id <- gff.mir.m@elementMetadata$ID[i]
        gtf.mir.m@elementMetadata$transcript_name <- gff.mir.m@elementMetadata$Name[i]
        gtf.mir.m@elementMetadata$exon_id[2] <- paste0(gff.mir.m@elementMetadata$ID[i],".e")
        gtf.mir.m@elementMetadata$gene_name <-
            gff.mir.p@elementMetadata$Name[gff.mir.p@elementMetadata$ID ==
                                           gff.mir.m@elementMetadata$Derives_from[i]]

        return(gtf.mir.m)}, mc.cores = nc)

    ## Concatenate gtf.mir.m 
    gtf.mir.m <- gtf.template
    slot(gtf.mir.m,"seqnames") <- do.call("c", mclapply(out, function(x) x@seqnames, mc.cores = nc))
    slot(gtf.mir.m,"strand") <- do.call("c", mclapply(out, function(x) x@strand, mc.cores = nc))
    slot(gtf.mir.m,"ranges") <- do.call("c", mclapply(out, function(x) x@ranges, mc.cores = nc))
    slot(gtf.mir.m,"elementMetadata")  <- do.call("rbind", mclapply(out, function(x) x@elementMetadata, mc.cores = nc))

    ## Modify transcript_name to avoid duplicates
    gtf.mir.m.t <- subset(gtf.mir.m, type == 'transcript')
    tid.dup <- gtf.mir.m.t@elementMetadata$transcript_name[
                 which(duplicated(gtf.mir.m.t@elementMetadata$transcript_name) == TRUE)]
    tidx.d <- which(gtf.mir.m@elementMetadata$transcript_name %in% tid.dup)

    suff <- sapply(tidx.d, function(i)
        sub(as.character(sub("-3p","",sub("-5p","", sub(
              "R","r",as.character(gtf.mir.m@elementMetadata$transcript_name[i]))))),
            "", gtf.mir.m@elementMetadata$gene_name[i]))
    gtf.mir.m@elementMetadata$transcript_name[tidx.d] <- paste0(
        gtf.mir.m@elementMetadata$transcript_name[tidx.d],suff)

    ## --------------------
    ## ---- 3. Re-format piRNA from RNACentral .gff annotations to .gtf

    ##  Get piRNA from RNAcentral database
    gff.pi <- gff.rcent[gff.rcent@elementMetadata$type.1 == "piRNA",]
    gff.pi <- gff.pi[as.character(gff.pi@seqnames) %in% c(as.character(c(1:19)),"X","Y","MT"),] 

    ## Remove pi annotations above 100 nt
    summary(gff.pi@ranges@width)
    gff.pi <- gff.pi[which(gff.pi@ranges@width < 100),]

    ## Find overlapping piRNA transcripts
    transcripts <- gff.pi[gff.pi@elementMetadata$type == "transcript",]
    ov <- IRanges::findOverlaps(transcripts,transcripts,type = "any", round(median(transcripts@ranges@width)/2))

    ## Select non-overlapping piRNA transcripts
    vnames <- transcripts@elementMetadata$ID
    g <- igraph::graph_from_data_frame(data.frame(
        from = vnames[ov@from], to = vnames[ov@to]),vertices = vnames)
    clu <- igraph::components(g)
    summary(clu$csize)
    clu.single <- names(clu$membership)[which(clu$membership %in% which(clu$csize==1))]

    ## Select one for overlapping piRNA annotations
    clu.primary <- mclapply(which(clu$csize!=1), function(i){

        ## Extract overlapping ids
        cl <-  names(clu$membership)[clu$membership == i]

        ## Find median transcript coordinates
        cl.df <- transcripts[which(transcripts@elementMetadata$ID %in% cl),]
        m.start <- median(cl.df@ranges@start)
        m.end <- median(cl.df@ranges@start + cl.df@ranges@width)

        primaryT <- which.min(abs(cl.df@ranges@start - m.start) +
                          abs(cl.df@ranges@start + cl.df@ranges@width - m.end))
        return(cl.df@elementMetadata$ID[primaryT])}, mc.cores = 30)

    gff.pi <- gff.pi[which(gff.pi@elementMetadata$ID %in%
                           c(unlist(clu.single),unlist(clu.primary)) |
                           as.character(gff.pi@elementMetadata$Parent) %in%
                           c(unlist(clu.single),unlist(clu.primary))),]

    ## Rename sequences
    gff.pi@seqnames <- droplevels(gff.pi@seqnames)

    ## Create an Ensemble like gtf empty templage
    gtf.template <- gtf.ens[1:3,]
    gtf.template@elementMetadata[1:3,-2] <- NA
    gtf.template@elementMetadata$gene_biotype <- rep("piRNA", 3)
    gtf.template@elementMetadata$transcript_biotype <- c(NA, rep("piRNA", 2))
    gtf.template@elementMetadata$source <- "RNACentral"

    ## piRNA .gtf generation loop
    idx <- which(gff.pi@elementMetadata$type == "transcript")
    out <- mclapply(idx, function(i) {

        if(i%%1000 == 0){print(i/tail(idx,1)*100)}
        
         ## Initialize .gtf template
         gtf.pi <- gtf.template

         ## Define annot coordinates
         gtf.pi@ranges[c(1:3)] <- gff.pi@ranges[i]
         gtf.pi@seqnames <- c(gff.pi@seqnames[i], gff.pi@seqnames[i], gff.pi@seqnames[i])
         gtf.pi@strand <- c(gff.pi@strand[i], gff.pi@strand[i], gff.pi@strand[i])

        ## Define metadata 
        gtf.pi@elementMetadata$gene_id[1:3] <- paste0(gff.pi@elementMetadata$ID[i],".g")
        gtf.pi@elementMetadata$transcript_id[2:3] <- gff.pi@elementMetadata$ID[i]
        gtf.pi@elementMetadata$exon_id[3] <- paste0(gff.pi@elementMetadata$ID[i],".e")
        gtf.pi@elementMetadata$gene_name[1:3] <- gtf.pi@elementMetadata$transcript_name[2:3] <-
            sub("_10090", "", gff.pi@elementMetadata$ID[i])

        return(gtf.pi)}, mc.cores = nc)

    ## Concatenate gtf.pi
    gtf.pi <- gtf.template
    tt <- mclapply(out, function(x) x@seqnames, mc.cores = nc); slot(gtf.pi,"seqnames") <- do.call("c", tt)
    tt <- mclapply(out, function(x) x@strand, mc.cores = nc); slot(gtf.pi,"strand") <- do.call("c", tt)
    tt <- mclapply(out, function(x) x@ranges, mc.cores = nc); slot(gtf.pi,"ranges") <- do.call("c", tt)
    tt <-  mclapply(out, function(x) x@elementMetadata, mc.cores = nc)
    slot(gtf.pi,"elementMetadata")  <- do.call("rbind",tt)

    ## --------------------
    ## ---- 4. Merge all .gtf annotations

    print("Exporting integrated annotations...")
    
    gtf.tmp <- gtf
    gtf.tmp@elementMetadata <- rbind(gtf.tmp@elementMetadata, gtf.pi@elementMetadata,
                                     gtf.mir.p@elementMetadata, gtf.mir.m@elementMetadata)
    gtf.tmp@seqnames <- c(gtf.tmp@seqnames, gtf.pi@seqnames, gtf.mir.p@seqnames, gtf.mir.m@seqnames)
    gtf.tmp@strand <- c(gtf.tmp@strand, gtf.pi@strand, gtf.mir.p@strand, gtf.mir.m@strand)
    gtf.tmp@ranges <- c(gtf.tmp@ranges, gtf.pi@ranges, gtf.mir.p@ranges, gtf.mir.m@ranges)
    gtf <- gtf.tmp
    
    ## Export annotations
    rtracklayer::export(gtf, paste0(db.dir,"/Mus_musculus.GRCm38.custom.gtf"))
    print("Mouse custom gtf successfully generated!")
}


## ARABIDOPSIS ANNOTATIONS EXTRACTION
## ----------------------------------------------------------------------------

if(annotArabidopsis){

    print("Importing annotations...")

    ## ---- Download annotations
    download.file(url = paste0("ftp://ftp.ensemblgenomes.org/pub/plants/current/gtf/",
                               "arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.50.gtf.gz"),
                  destfile = paste0(db.dir,"ensembl_Arabidopsis_thaliana.TAIR10.50.gtf.gz"))
    download.file(url = "ftp://mirbase.org/pub/mirbase/CURRENT/genomes/ath.gff3",
                  destfile = paste0(db.dir,"mirbase_ath.gff3"))
    download.file(url = paste0("ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/16.0/",
                               "genome_coordinates/gff3/arabidopsis_thaliana.TAIR10.gff3.gz"),
                  destfile = paste0(db.dir,"rnacentral_arabidopsis_thaliana.TAIR10.gff3.gz"))
                  
    ## ---- Import raw annotations
    print("Importing annotations...")
    gtf.ens <- import(paste0(db.dir,"ensembl_Arabidopsis_thaliana.TAIR10.50.gtf.gz"))
    gff.mir <- import(paste0(db.dir,"mirbase_ath.gff3"))
    gff.rcent <- import(paste0(db.dir,"rnacentral_arabidopsis_thaliana.TAIR10.gff3.gz"))

    ## Add transcript_name missing column (for Ensembl, transcript_name = transcript_id)
    gtf.ens@elementMetadata$transcript_name <- gtf.ens@elementMetadata$transcript_id
    
    ## Exclude miRNA annotations
    gtf <- gtf.ens[gtf.ens@elementMetadata$gene_biotype != "miRNA",]
    
    ## ---- 2. Re-format pre-miRNA miRbase gff annotations to .gtf
    print("Integrating miRNA from miRbase...")
    
    ## Relabel sequence names
    gff.mir@seqnames@values <- mapvalues(gff.mir@seqnames@values, from = levels(gff.mir@seqnames),
                                         to = c(as.character(c(1:length(levels(gff.mir@seqnames))))))

    gff.mir.p <- gff.mir[gff.mir@elementMetadata$type == "miRNA_primary_transcript"]
    gff.mir.m <- gff.mir[gff.mir@elementMetadata$type == "miRNA"]

    ## --------------------
    ## pre-miRNA annotations

    ## Create an Ensemble like gtf empty templage
    gtf.template <- gtf.ens[1,]
    gtf.template@elementMetadata[,-2] <- NA
    gtf.template@elementMetadata$gene_biotype <- "pre_miRNA"
    gtf.template@elementMetadata$transcript_biotype <- NA
    gtf.template@elementMetadata$source <- "mirbase"
    
    ## pre-miRNA .gtf generation loop
    out <- mclapply(c(1:length(gff.mir.p@seqnames)), function(i) {

        ## Initialize .gtf template
        gtf.mir.p <- gtf.template

        ## Define annot coordinates
        gtf.mir.p@ranges <- gff.mir.p@ranges[i]
        gtf.mir.p@seqnames <- gff.mir.p@seqnames[i]
        gtf.mir.p@strand <- gff.mir.p@strand[i]

        ## Define metadata 
        gtf.mir.p@elementMetadata$gene_id <- gff.mir.p@elementMetadata$ID[i]
        gtf.mir.p@elementMetadata$gene_name <- gff.mir.p@elementMetadata$Name[i]

        return(gtf.mir.p)}, mc.cores = nc)

    ## Concatenate gtf.mir.p 
    gtf.mir.p <- gtf.template
    slot(gtf.mir.p,"seqnames") <- do.call("c", mclapply(out, function(x) x@seqnames, mc.cores = nc))
    slot(gtf.mir.p,"strand") <- do.call("c", mclapply(out, function(x) x@strand, mc.cores = nc))
    slot(gtf.mir.p,"ranges") <- do.call("c", mclapply(out, function(x) x@ranges, mc.cores = nc))
    slot(gtf.mir.p,"elementMetadata")  <- do.call("rbind", mclapply(out, function(x) x@elementMetadata, mc.cores = nc))

    ## --------------------
    ## Mature miRNA annotations
        
    ## Create an Ensemble like gtf empty templage
    gtf.template <- gtf.ens[1:2,]
    gtf.template@elementMetadata[,] <- NA
    gtf.template@elementMetadata$gene_biotype <- "pre_miRNA"
    gtf.template@elementMetadata$transcript_biotype <- "miRNA"
    gtf.template@elementMetadata$source <- "mirbase"
    gtf.template@elementMetadata$type <- c("transcript","exon")

    ## mature-miRNA .gtf generation loop
    out <- mclapply(c(1:length(gff.mir.m@seqnames)), function(i) {

        ## Initialize .gtf template
        gtf.mir.m <- gtf.template

        ## Define annot coordinates
        gtf.mir.m@ranges[c(1,2)] <- gff.mir.m@ranges[i]
        gtf.mir.m@seqnames <- c(gff.mir.m@seqnames[i], gff.mir.m@seqnames[i])
        gtf.mir.m@strand <- c(gff.mir.m@strand[i], gff.mir.m@strand[i])

        ## Define metadata
        gtf.mir.m@elementMetadata$gene_id <- gff.mir.m@elementMetadata$Derives_from[i]
        gtf.mir.m@elementMetadata$transcript_id <- gff.mir.m@elementMetadata$ID[i]
        gtf.mir.m@elementMetadata$transcript_name <- gff.mir.m@elementMetadata$Name[i]
        gtf.mir.m@elementMetadata$exon_id[2] <- paste0(gff.mir.m@elementMetadata$ID[i],".e")
        gtf.mir.m@elementMetadata$gene_name <-
            gff.mir.p@elementMetadata$Name[gff.mir.p@elementMetadata$ID ==
                                           gff.mir.m@elementMetadata$Derives_from[i]]

        return(gtf.mir.m)}, mc.cores = nc)

    ## Concatenate gtf.mir.m 
    gtf.mir.m <- gtf.template
    slot(gtf.mir.m,"seqnames") <- do.call("c", mclapply(out, function(x) x@seqnames, mc.cores = nc))
    slot(gtf.mir.m,"strand") <- do.call("c", mclapply(out, function(x) x@strand, mc.cores = nc))
    slot(gtf.mir.m,"ranges") <- do.call("c", mclapply(out, function(x) x@ranges, mc.cores = nc))
    slot(gtf.mir.m,"elementMetadata")  <- do.call("rbind", mclapply(out, function(x)
        x@elementMetadata, mc.cores = nc))

    
    ## --------------------
    ## ---- 3. Re-format siRNA from RNACentral .gff annotations to .gtf

    ##  Get siRNA from RNAcentral database
    gff.si <- gff.rcent[gff.rcent@elementMetadata$type.1 == "siRNA",]
    levels(gtf.ens@seqnames) == levels(gff.si@seqnames)

    ## Remove siRNA annotations above 100 nt
    summary(gff.si@ranges@width)
    gff.si <- gff.si[which(gff.si@ranges@width < 100),]

    ## Find overlapping siRNA transcripts, first remove siRNA overlapping with miRNA (
    ## redundant annotation of microRNA as a siRNA)
    transcripts <- gff.si[gff.si@elementMetadata$type == "transcript",]
    ov.mi <- IRanges::findOverlaps(transcripts, gff.rcent[gff.rcent@elementMetadata$type.1 == "miRNA" &
                                                          gff.rcent@elementMetadata$type == 'transcript',],
                                   type = "any", round(median(transcripts@ranges@width)/2))
    transcripts <- transcripts[-ov.mi@from,]
    
    ## Find overlapping siRNA transcripts
    ov <- IRanges::findOverlaps(transcripts,transcripts,type = "any", round(median(transcripts@ranges@width)/2))

    ## Select non-overlapping piRNA transcripts
    vnames <- transcripts@elementMetadata$ID
    g <- igraph::graph_from_data_frame(data.frame(
        from = vnames[ov@from], to = vnames[ov@to]), vertices = vnames)
    clu <- igraph::components(g)
    summary(clu$csize)
    clu.single <- names(clu$membership)[which(clu$membership %in% which(clu$csize==1))]

    ## Select one for overlapping siRNA annotations
    clu.primary <- mclapply(which(clu$csize!=1), function(i){

        ## Extract overlapping ids
        cl <-  names(clu$membership)[clu$membership == i]

        ## Find median transcript coordinates
        cl.df <- transcripts[which(transcripts@elementMetadata$ID %in% cl),]
        m.start <- median(cl.df@ranges@start)
        m.end <- median(cl.df@ranges@start + cl.df@ranges@width)

        primaryT <- which.min(abs(cl.df@ranges@start - m.start) +
                          abs(cl.df@ranges@start + cl.df@ranges@width - m.end))
        return(cl.df@elementMetadata$ID[primaryT])}, mc.cores = 30)

    gff.si <- gff.si[which(gff.si@elementMetadata$ID %in%
                           c(unlist(clu.single),unlist(clu.primary)) |
                           as.character(gff.si@elementMetadata$Parent) %in%
                           c(unlist(clu.single),unlist(clu.primary))),]

    ## Rename sequences
    gff.si@seqnames <- droplevels(gff.si@seqnames)

    ## Create an Ensemble like gtf empty templage
    gtf.template <- gtf.ens[1:3,]
    gtf.template@elementMetadata[1:3,-2] <- NA
    gtf.template@elementMetadata$gene_biotype <- rep("siRNA", 3)
    gtf.template@elementMetadata$transcript_biotype <- c(NA, rep("siRNA", 2))
    gtf.template@elementMetadata$source <- "RNACentral"

    ## siRNA .gtf generation loop
    idx <- which(gff.si@elementMetadata$type == "transcript")
    out <- mclapply(idx, function(i) {

        if(i%%1000 == 0){print(i/tail(idx,1)*100)}
        
         ## Initialize .gtf template
         gtf.si <- gtf.template

         ## Define annot coordinates
         gtf.si@ranges[c(1:3)] <- gff.si@ranges[i]
         gtf.si@seqnames <- c(gff.si@seqnames[i], gff.si@seqnames[i], gff.si@seqnames[i])
         gtf.si@strand <- c(gff.si@strand[i], gff.si@strand[i], gff.si@strand[i])

        ## Define metadata 
        gtf.si@elementMetadata$gene_id[1:3] <- paste0(gff.si@elementMetadata$Name[i],".g")
        gtf.si@elementMetadata$transcript_id[2:3] <- gff.si@elementMetadata$Name[i]
        gtf.si@elementMetadata$exon_id[3] <- paste0(gff.si@elementMetadata$Name[i],".e")

        return(gtf.si)}, mc.cores = nc)

    ## Concatenate gtf.si
    gtf.si <- gtf.template
    slot(gtf.si,"seqnames") <- do.call("c", mclapply(out, function(x) x@seqnames, mc.cores = nc))
    slot(gtf.si,"strand") <- do.call("c", mclapply(out, function(x) x@strand, mc.cores = nc))
    slot(gtf.si,"ranges") <- do.call("c", mclapply(out, function(x) x@ranges, mc.cores = nc))
    slot(gtf.si,"elementMetadata")  <- do.call("rbind", mclapply(out, function(x) x@elementMetadata, mc.cores = nc))
    

    ## --------------------
    ## ---- 4. Merge all .gtf annotations

    print("Exporting integrated annotations...")
    
    gtf.tmp <- gtf
    gtf.tmp@elementMetadata <- rbind(gtf.tmp@elementMetadata, gtf.si@elementMetadata,
                                     gtf.mir.p@elementMetadata, gtf.mir.m@elementMetadata)
    gtf.tmp@seqnames <- c(gtf.tmp@seqnames, gtf.si@seqnames, gtf.mir.p@seqnames, gtf.mir.m@seqnames)
    gtf.tmp@strand <- c(gtf.tmp@strand, gtf.si@strand, gtf.mir.p@strand, gtf.mir.m@strand)
    gtf.tmp@ranges <- c(gtf.tmp@ranges, gtf.si@ranges, gtf.mir.p@ranges, gtf.mir.m@ranges)
    gtf <- gtf.tmp
    
    ## Export annotations
    rtracklayer::export(gtf, paste0(db.dir,"/Arabidopsis_thaliana.TAIR10.custom.gtf"))
    print("Arabidopsis custom gtf successfully generated!")
}


## NEMATODE ANNOTATIONS EXTRACTION
## ----------------------------------------------------------------------------

if(annotNematode){
    
    ## ## ---- Download annotations
    ## download.file(url = paste0("ftp://ftp.ensembl.org/pub/current_gtf/caenorhabditis_elegans",
    ##                            "/Caenorhabditis_elegans.WBcel235.103.gtf.gz"),
    ##               destfile = paste0(db.dir,"/ensembl_Caenorhabditis_elegans.WBcel235.103.gtf.gz"))
    
    ## download.file(url = paste0("ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/",
    ##                            "current/species/caenorhabditis_elegans/PRJNA13758/",
    ##                            "caenorhabditis_elegans.PRJNA13758.WBPS15.canonical_geneset.gtf.gz"),
    ##               destfile = paste0(db.dir,"wormbase_",
    ##                                 "caenorhabditis_elegans.PRJNA13758.WBPS15.canonical_geneset.gtf.gz"))

    ## download.file(url = "ftp://mirbase.org/pub/mirbase/CURRENT/genomes/cel.gff3",
    ##               destfile = paste0(db.dir,'mirbase_cel.gff3'))

    ## download.file(url = paste0("ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/16.0/",
    ##                            "genome_coordinates/gff3/caenorhabditis_elegans.WBcel235.gff3.gz"),
    ##               destfile = paste0(db.dir,"rnacentral_caenorhabditis_elegans.WBcel235.gff3.gz"))

    ## ---- Import raw annotations
    print("Importing annotations...")
    
    gtf.ens <- import(paste0(db.dir,"ensembl_Caenorhabditis_elegans.WBcel235.103.gtf.gz"))
    gtf.worm <- import(paste0(db.dir,"wormbase_caenorhabditis_elegans.PRJNA13758.WBPS15.canonical_geneset.gtf.gz"))
    gff.rcent <- import(paste0(db.dir,"/rnacentral_caenorhabditis_elegans.WBcel235.gff3.gz"))
    gff.mir <- import(paste0(db.dir, 'mirbase_cel.gff3'))

    ## ---- Explore databases
    with(subset(gtf.ens@elementMetadata, type == "gene"), table(gene_biotype))
    with(subset(gtf.worm@elementMetadata, type == "gene"), table(gene_biotype))
    with(gff.mir@elementMetadata, table(type))
    with(subset(gff.rcent@elementMetadata, type == "transcript"), table(type.1))
    with(subset(gtf.ens@elementMetadata, type == "gene"), table(duplicated(gene_name)))

    ## Add transcript_name missing column (for Ensembl, set transcript_name = transcript_id)
    gtf.ens@elementMetadata$transcript_name <- gtf.ens@elementMetadata$transcript_id
    
    ## Exclude miRNA annotations from Ensembl
    gtf <- gtf.ens[gtf.ens@elementMetadata$gene_biotype != "miRNA",]
    
    ## ---- 2. Re-format pre-miRNA miRbase gff annotations to .gtf

    print("Integrating miRNA from miRbase...")
    
    ## Relabel sequence names
    gff.mir@seqnames@values <- mapvalues(gff.mir@seqnames@values, from = levels(gff.mir@seqnames),
                                         to =  sub('chr', '', levels(gff.mir@seqnames)))

    gff.mir.p <- gff.mir[gff.mir@elementMetadata$type == "miRNA_primary_transcript"]
    gff.mir.m <- gff.mir[gff.mir@elementMetadata$type == "miRNA"]

    ## --------------------
    ## pre-miRNA annotations

    ## Create an Ensemble like gtf empty templage
    gtf.template <- gtf.ens[1,]
    gtf.template@elementMetadata[,-2] <- NA
    gtf.template@elementMetadata$gene_biotype <- "pre_miRNA"
    gtf.template@elementMetadata$transcript_biotype <- NA
    gtf.template@elementMetadata$source <- "mirbase"
    
    ## pre-miRNA .gtf generation loop
    out <- mclapply(c(1:length(gff.mir.p@seqnames)), function(i) {

        ## Initialize .gtf template
        gtf.mir.p <- gtf.template

        ## Define annot coordinates
        gtf.mir.p@ranges <- gff.mir.p@ranges[i]
        gtf.mir.p@seqnames <- gff.mir.p@seqnames[i]
        gtf.mir.p@strand <- gff.mir.p@strand[i]

        ## Define metadata 
        gtf.mir.p@elementMetadata$gene_id <- gff.mir.p@elementMetadata$ID[i]
        gtf.mir.p@elementMetadata$gene_name <- gff.mir.p@elementMetadata$Name[i]

        return(gtf.mir.p)}, mc.cores = nc)

    ## Concatenate gtf.mir.p 
    gtf.mir.p <- gtf.template
    slot(gtf.mir.p,"seqnames") <- do.call("c", mclapply(out, function(x) x@seqnames, mc.cores = nc))
    slot(gtf.mir.p,"strand") <- do.call("c", mclapply(out, function(x) x@strand, mc.cores = nc))
    slot(gtf.mir.p,"ranges") <- do.call("c", mclapply(out, function(x) x@ranges, mc.cores = nc))
    slot(gtf.mir.p,"elementMetadata")  <- do.call("rbind", mclapply(out, function(x) x@elementMetadata, mc.cores = nc))

    ## --------------------
    ## Mature miRNA annotations
        
    ## Create an Ensemble like gtf empty templage
    gtf.template <- gtf.ens[1:2,]
    gtf.template@elementMetadata[,] <- NA
    gtf.template@elementMetadata$gene_biotype <- "pre_miRNA"
    gtf.template@elementMetadata$transcript_biotype <- "miRNA"
    gtf.template@elementMetadata$source <- "mirbase"
    gtf.template@elementMetadata$type <- c("transcript","exon")

    ## mature-miRNA .gtf generation loop
    out <- mclapply(c(1:length(gff.mir.m@seqnames)), function(i) {

        ## Initialize .gtf template
        gtf.mir.m <- gtf.template

        ## Define annot coordinates
        gtf.mir.m@ranges[c(1,2)] <- gff.mir.m@ranges[i]
        gtf.mir.m@seqnames <- c(gff.mir.m@seqnames[i], gff.mir.m@seqnames[i])
        gtf.mir.m@strand <- c(gff.mir.m@strand[i], gff.mir.m@strand[i])

        ## Define metadata
        gtf.mir.m@elementMetadata$gene_id <- gff.mir.m@elementMetadata$Derives_from[i]
        gtf.mir.m@elementMetadata$transcript_id <- gff.mir.m@elementMetadata$ID[i]
        gtf.mir.m@elementMetadata$transcript_name <- gff.mir.m@elementMetadata$Name[i]
        gtf.mir.m@elementMetadata$exon_id[2] <- paste0(gff.mir.m@elementMetadata$ID[i],".e")
        gtf.mir.m@elementMetadata$gene_name <-
            gff.mir.p@elementMetadata$Name[gff.mir.p@elementMetadata$ID ==
                                           gff.mir.m@elementMetadata$Derives_from[i]]

        return(gtf.mir.m)}, mc.cores = nc)

    ## Concatenate gtf.mir.m 
    gtf.mir.m <- gtf.template
    slot(gtf.mir.m,"seqnames") <- do.call("c", mclapply(out, function(x) x@seqnames, mc.cores = nc))
    slot(gtf.mir.m,"strand") <- do.call("c", mclapply(out, function(x) x@strand, mc.cores = nc))
    slot(gtf.mir.m,"ranges") <- do.call("c", mclapply(out, function(x) x@ranges, mc.cores = nc))
    slot(gtf.mir.m,"elementMetadata")  <- do.call("rbind", mclapply(out, function(x)
        x@elementMetadata, mc.cores = nc))


    ## Modify transcript_name to avoid duplicates
    gtf.mir.m.t <- subset(gtf.mir.m, type == 'transcript')
    tid.dup <- gtf.mir.m.t@elementMetadata$transcript_name[
                 which(duplicated(gtf.mir.m.t@elementMetadata$transcript_name) == TRUE)]
    tidx.d <- which(gtf.mir.m@elementMetadata$transcript_name %in% tid.dup)

    suff <- sapply(tidx.d, function(i)
        sub(as.character(sub("-3p","",sub("-5p","", sub(
              "R","r",as.character(gtf.mir.m@elementMetadata$transcript_name[i]))))),
            "", gtf.mir.m@elementMetadata$gene_name[i]))
    gtf.mir.m@elementMetadata$transcript_name[tidx.d] <- paste0(
        gtf.mir.m@elementMetadata$transcript_name[tidx.d],suff)
        

    ## --------------------
    ## ---- 4. Merge all .gtf annotations

    print("Exporting integrated annotations...")
    
    gtf.tmp <- gtf
    gtf.tmp@elementMetadata <-rbind(gtf.tmp@elementMetadata,gtf.mir.p@elementMetadata,gtf.mir.m@elementMetadata)
    gtf.tmp@seqnames <- c(gtf.tmp@seqnames, gtf.mir.p@seqnames, gtf.mir.m@seqnames)
    gtf.tmp@strand <- c(gtf.tmp@strand, gtf.mir.p@strand, gtf.mir.m@strand)
    gtf.tmp@ranges <- c(gtf.tmp@ranges, gtf.mir.p@ranges, gtf.mir.m@ranges)
    gtf <- gtf.tmp

    with(subset(gtf@elementMetadata, type == "gene"), table(gene_biotype))
    with(subset(gtf@elementMetadata, type == "gene"), table(duplicated(gene_name)))

    with(subset(gtf@elementMetadata, type == "transcript"), table(transcript_biotype))
    with(subset(gtf@elementMetadata, type == "transcript"), table(duplicated(transcript_name)))
    
    ## Export annotations
    rtracklayer::export(gtf, paste0(db.dir,"/Caenorhabditis_elegans.WBcel235.custom.gtf"))
    print("Nematode custom gtf successfully generated!")
}


