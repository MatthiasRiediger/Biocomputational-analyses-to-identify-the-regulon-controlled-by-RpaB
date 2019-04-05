# Required packages
require(plyr)
require(reshape2)
require(seqinr)
require(Biostrings) # Bioconductor package
require(TFBSTools) # Bioconductor package
require(ggplot2)


# Input 1: Import promoter database (from 400nt upstream to 200nt downstream of TSS). Information derived from Kopf et al., 2014	
promoter <- readDNAStringSet("promoter_database.txt")	
  promoter_length <- nchar(promoter[1])-201

  ACGT <- letterFrequency(promoter, letters="ACGT", OR=0)							
    sumACGT <- sum(ACGT)
    sumA <- sum(ACGT[,1]) / sumACGT		
    sumC <- sum(ACGT[,2]) / sumACGT	
    sumG <- sum(ACGT[,3]) / sumACGT
    sumT <- sum(ACGT[,4]) / sumACGT	


# Input 2: Import annotation file including comparative dRNA-seq data. Information derived from Kopf et al., 2014
annotation <- read.csv(file="Annotationfile.csv", header = TRUE, sep = ",", quote = "\"", dec = ".", fill = TRUE)


# Input 3: Import TF binding motif (HLR1) from published binding sites. Information derived from Riediger et al., 2018
input_motif <- read.csv("HLR1.txt") 
  motif_name <- colnames(input_motif)	



# Function: perform motif search using SearchSeq 
motif_search <- function(motif_name, PFM, sumA, sumC, sumG, sumT, promoter, promoter_length)	{	 
  # Transform PFM to PWM 	
  PWM <- toPWM(PFM, type="log2probratio", pseudocounts=0.8, bg=c(A=sumA, C=sumC, G=sumG, T=sumT))

  # Apply SearchSeq function to search for input_motif in promoter database (TFBSTools package)
  raw_output <- searchSeq(PWM, promoter, strand="*", min.score="70%")
  raw_output <- transform(as(raw_output, "data.frame"), input = motif_name)
  raw_output <- raw_output[order(-raw_output$relScore),]

  # Limit results to max 5 motifs per promoter 
  raw_output <- 	transform(raw_output, promoterspecific_motifrank=ave(1:nrow(raw_output), seqnames,FUN=function(x) order(relScore[x],decreasing=TRUE)))
  raw_output <- droplevels(subset(raw_output, promoterspecific_motifrank <= 5))

  # Format raw_output 
  raw_output$start <- raw_output$start-promoter_length
  raw_output$end <- raw_output$end-promoter_length
  raw_output <- cbind(raw_output[,c("seqnames", "start", "end", "absScore", "relScore", "strand","siteSeqs", "input", "promoterspecific_motifrank"),drop=FALSE])
    names(raw_output)[names(raw_output) == "siteSeqs"] <- 'motifSeq'				
    names(raw_output)[names(raw_output) == "seqnames"] <- 'TU.ID'	
  raw_output <<- raw_output
}



# Function: calculate the average relative expression per position of the motif
AvRelExpr <- function(output, input_motif){
  
  half_motif_length <- nchar(as.character(input_motif[1,1]))/2
  
  sub_output <- droplevels(subset(output, relScore >= 0.9))
  sub_output <- sub_output[order(-sub_output$relScore),]
  
  sub_av_output <- c()	
  for (i in 1:length(levels(sub_output$input))){				
    x.sub <- c()
    y.sub <- sub_output[which(sub_output$input == levels(sub_output$input)[i]),]
    y.sub$start <- as.factor(y.sub$start)
    y.sub[] <- lapply(y.sub, function(x) if(is.factor(x)) factor(x) else x)
    for(j in 1:length(levels(y.sub$start))){										
      z.sub <- which(y.sub$start == levels(y.sub$start)[j]) 
      z.sub <- y.sub[z.sub,]
      z.sub[] <- lapply(z.sub, function(x) if(is.factor(x)) factor(x) else x)
      z.sub <- transform(z.sub, 'av 15C/Median' = mean((log2(z.sub$X15C/z.sub$Median)+1))-1)
      z.sub <- transform(z.sub, 'av CO2/Median' = mean((log2(z.sub$X.CO2/z.sub$Median)+1))-1)
      z.sub <- transform(z.sub, 'av 42C/Median' = mean((log2(z.sub$X42C/z.sub$Median)+1))-1)
      z.sub <- transform(z.sub, 'av Darkness/Median' = mean((log2(z.sub$Darkness/z.sub$Median)+1))-1)
      z.sub <- transform(z.sub, 'av Fe/Median' = mean((log2(z.sub$X.Fe/z.sub$Median)+1))-1)
      z.sub <- transform(z.sub, 'av HL/Median' = mean((log2(z.sub$Hl/z.sub$Median)+1))-1)
      z.sub <- transform(z.sub, 'av N/Median' = mean((log2(z.sub$X.N/z.sub$Median)+1))-1)
      z.sub <- transform(z.sub, 'av P/Median' = mean((log2(z.sub$X.P/z.sub$Median)+1))-1)
      z.sub <- transform(z.sub, 'av ExpP/Median' = mean((log2(z.sub$Exp..p./z.sub$Median)+1))-1)
      z.sub <- transform(z.sub, 'av StatP/Median' = mean((log2(z.sub$Stat..p./z.sub$Median)+1))-1)
      x.sub <- rbind(x.sub, z.sub)}
    sub_av_output <- rbind(sub_av_output, x.sub)
  } 
  

  sub_av_output[] <- lapply(sub_av_output, function(x) if(is.factor(x)) factor(x) else x)	
  sub_av_output <- cbind(sub_av_output[,c("input", "start"),drop=FALSE], sub_av_output[,c('av.15C.Median', 
                              'av.CO2.Median', 'av.42C.Median' , 'av.Darkness.Median' , 'av.Fe.Median', 'av.HL.Median', 
                              'av.N.Median', 'av.P.Median', 'av.ExpP.Median', 'av.StatP.Median'),drop=FALSE])


  sub_av_output <- sub_av_output[!duplicated(sub_av_output[,1:2]),]
  sub_av_output <- melt(sub_av_output, id.vars=c("input", "start"))
  names(sub_av_output)[3:4] <- c("condition", "av.rel.Expression")
  sub_av_output$av.rel.Expression[is.infinite(sub_av_output$av.rel.Expression)] <- 0
  sub_av_output$av.rel.Expression[is.na(sub_av_output$av.rel.Expression)] <- 0
  sub_av_output$av.rel.Expression[is.nan(sub_av_output$av.rel.Expression)] <- 0
  sub_av_output <- cbind(sub_av_output[,2:4], data.frame(do.call(rbind, strsplit(as.vector(sub_av_output$input), split = "_"))))
  names(sub_av_output)[4:5] <- c("input", "number")
    AvExprPerPos <<- sub_av_output


  # Figure S2
p1<-ggplot(sub_av_output, aes(x=as.numeric(levels(start)[start])+ half_motif_length, y=av.rel.Expression, color=input)) + 
    geom_smooth(mapping = NULL, data = NULL, stat = "smooth", position = "identity", method = "auto", se = TRUE, n = 10,level = 0.25, na.rm=TRUE) +
    scale_x_continuous(name = "Distance to TSS", expand=c(0,0), breaks=c(200,100,50,0,-50,-100,-200,-400)) +
    scale_y_continuous(name="Mean log2[fc]",expand=c(0,0))	+	
    coord_cartesian(ylim=c(2,-2))+
    theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(),
          legend.background = element_rect(fill="white", size=0.05, linetype="solid", colour ="black"), 
          legend.title = element_blank(), legend.text=element_text(size=8), axis.line.x = element_line(color="black", size = .1), 
          axis.line.y = element_line(color="black", size = .1), axis.text=element_text(size=8),
          axis.title=element_text(size=10)) +
    geom_vline(aes(xintercept=0), colour="black", linetype="solid") +
    geom_hline(aes(yintercept=0), colour="black", linetype="solid") +
    facet_wrap( ~ condition, nrow=5)
  p1 
  ggsave("AverageExpressionPerPosition.pdf", width=12, height=16, dpi=1200)

} 



# Function: define enriched motif areas by comparison to background using a density threshold
EnrichedAreas <- function(output, input_motif) {
  half_motif_length <- nchar(as.character(input_motif[1,1]))/2 
  width <- 5
  
  background <- droplevels(subset(output, input != motif_name))
  density_threshold <- quantile(density(background$start, bw=width)$y[0:length(density(background$start, bw=width)$y)], c(.99))[1]
  
  sub_output <- droplevels(subset(output, relScore >= 0.9))
  act_motif <- droplevels(subset(sub_output, input == motif_name & Hl/Median < 1))
  repr_motif <- droplevels(subset(sub_output, input == motif_name & Hl/Median >= 1))
  act_background <- droplevels(subset(sub_output, input != motif_name & Hl/Median < 1))
  repr_background <- droplevels(subset(sub_output, input != motif_name & Hl/Median >= 1))
  
  act_motif_density <- as.data.frame(cbind(round(density(act_motif$start, bw= width)$x,0), density(act_motif$start, bw= width)$y))
  names(act_motif_density)[1:2] <- c("start", "act_motif_density")
  repr_motif_density <- as.data.frame(cbind(round(density(repr_motif$start, bw= width)$x,0), density(repr_motif$start, bw= width)$y))
  names(repr_motif_density)[1:2] <- c("start", "repr_motif_density")
  act_background_density <- as.data.frame(cbind(round(density(act_background$start, bw= width)$x,0), density(act_background$start, bw= width)$y))
  names(act_background_density)[1:2] <- c("start", "act_background_density")
  repr_background_density <- as.data.frame(cbind(round(density(repr_background$start, bw= width)$x,0), density(repr_background$start, bw= width)$y))
  names(repr_background_density)[1:2] <- c("start", "repr_background_density")
  
  act_density_comparison <- merge(act_motif_density, act_background_density, by="start")
  act_density_comparison <- transform(act_density_comparison, center = start + half_motif_length)
  act_density_comparison <- transform(act_density_comparison, threshold = density_threshold )
  
  repr_density_comparison <- merge(repr_motif_density, repr_background_density, by="start")
  repr_density_comparison <- transform(repr_density_comparison, center = start + half_motif_length)
  repr_density_comparison <- transform(repr_density_comparison, threshold = density_threshold )
  
  act_enriched_area <<- droplevels(subset(act_density_comparison, act_motif_density - threshold >= 0 & act_motif_density - act_background_density >= 0))
  repr_enriched_area <<- droplevels(subset(repr_density_comparison, repr_motif_density - threshold >= 0 & repr_motif_density - repr_background_density >= 0))
  
}



# Calculate probabilites for each motif from promoter position (enriched area) and relative Expression
RankProbability <- function(act_center_max, act_center_min, repr_center_max, repr_center_min, output, input_motif, motif_name, n_background){
  
  half_motif_length <- nchar(as.character(input_motif[1,1]))/2
  
  output <- output[order(-output$relScore),]
  
  act <- droplevels(subset(output, start <= act_center_min-half_motif_length & start >= act_center_max-half_motif_length & relScore > 0.8)) #1st subset: activated location only
  repr <- droplevels(subset(output, start <= repr_center_min-half_motif_length & start >= repr_center_max-half_motif_length & relScore > 0.8)) #1st subset: repressed location only
  
  
  act <- transform(act, score_rank = rank(-relScore, na.last = "keep", ties.method = "first"))
  act <- transform(act, expr_rank_repr = rank(-Hl/Median, na.last = "keep", ties.method = "first"))
  act <- transform(act, expr_rank_act = rank(Hl/Median, na.last = "keep", ties.method = "first"))
  act <- transform(act, combined_rank_total = rank((score_rank + expr_rank_act), na.last = "keep", ties.method = "first"))
  act <- transform(act, combined_rank_motif = ave(combined_rank_total, input, FUN= function(combined_rank_total) rank(combined_rank_total, na.last = "keep", ties.method = "first")))
  act <- transform(act, probability = 1-((combined_rank_total - combined_rank_motif)/n_background/combined_rank_motif))
  bg_act <- droplevels(subset(act, input != motif_name))
  act <- droplevels(subset(act, input == motif_name))
  act <- act[order(act$combined_rank_motif),]
  
  repr <- transform(repr, score_rank = rank(-relScore, na.last = "keep", ties.method = "first"))
  repr <- transform(repr, expr_rank_repr = rank(-Hl/Median, na.last = "keep", ties.method = "first"))
  repr <- transform(repr, expr_rank_act = rank(Hl/Median, na.last = "keep", ties.method = "first"))
  repr <- transform(repr, combined_rank_total = rank((score_rank + expr_rank_repr), na.last = "keep", ties.method = "first"))
  repr <- transform(repr, combined_rank_motif = ave(combined_rank_total, input, FUN= function(combined_rank_total) rank(combined_rank_total, na.last = "keep", ties.method = "first")))
  repr <- transform(repr, probability = 1-((combined_rank_total - combined_rank_motif)/n_background/combined_rank_motif))
  bg_repr <- droplevels(subset(repr, input != motif_name))
  repr <- droplevels(subset(repr, input == motif_name))
  repr <- repr[order(repr$combined_rank_motif),]
  
  
  final_output <<- rbind(act, repr)
}



# Transform input_motif to PFM 	
PFM <- consensusMatrix(as.character(input_motif[,1]))	
  PFM <- PFMatrix(ID="", name="", matrixClass="", strand="+", bg=c(A=sumA, C=sumC, G=sumG, T=sumT), 
                tags=list(family="", species="", tax_group="",medline="", type="",ACC="", pazar_tf_id="",
                TFBSshape_ID="", TFencyclopedia_ID=""), profileMatrix=PFM)	



# Motif search 
motif_search(motif_name, PFM, sumA, sumC, sumG, sumT, promoter, promoter_length) # output data.frame: raw_output
  motif_output <- raw_output


  
# Repeat motif search 100 times with shuffled input motif to generate background model
n_background <- 100 # Select no of motif shuffling


background_output <- c()
for(i in 1:n_background){
  # Shuffle input_motif --> Create background_motif
  PFM <- permuteMatrix(PFM, type="intra")						
  
  background_name <- paste0("backgroundmotif_", i)
  
  # Motif search
  motif_search(background_name, PFM, sumA, sumC, sumG, sumT, promoter, promoter_length)
  
  background_output <- rbind(background_output, raw_output)
} # output data.frame: background_output



# Annotation file merge
output <- rbind(motif_output, background_output)
output <- merge(output, annotation, by="TU.ID", all.y = FALSE)



#Calculate average relative expression per position of the motif
AvRelExpr(output, input_motif) # output data.frame: AvExprPerPos + Figure S2



# Define enriched motif areas by comparison to background using a density threshold
    # Selected condition: Hl (from 10 conditions: X15C, X42C, X.CO2, Darkness, X.Fe,	Hl = Highlight,	X.N, X.P, Exp..p., Stat..p.)
EnrichedAreas(output, input_motif) # output data.frame: act_enriched_area [=Highlight < 1], repr_enriched_area [=Highlight >= 1]



# Calculate probabilites for each motif from promoter position and relative Expression
act_center_max <- -66 # select upper end for activating area
act_center_min <- -45 # select lower end for activating area
repr_center_max <- -38 # select upper end for repressing area
repr_center_min <- +23 # select lower end for repressing area


RankProbability(act_center_max, act_center_min, repr_center_max, repr_center_min, output, input_motif, motif_name, n_background) # output data.frame: final_output
 


