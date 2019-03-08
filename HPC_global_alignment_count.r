## Author: Qi Yu
## Contact: dryuqi@gmail.com

library(seqinr)
library(Biostrings)
library(MASS)
library(plyr)
library(ggplot2)

args<-(commandArgs(TRUE))



a1<-read.alignment("HPC_4N_B6XCAST_Alu_new_10k.txt.fa_200",format="fasta")
#head(a1)

printPairwiseAlignment <- function(alignment, chunksize)
{
  seq1aln <- pattern(alignment) # Get the alignment for the firstsequence
  seq2aln <- subject(alignment) # Get the alignment for the secondsequence
  alnlen  <- nchar(seq1aln)     # Find the number of columns in thealignment
  starts  <- seq(1, alnlen, by=chunksize)
  n       <- length(starts)     
  seq1alnresidues <- 0
  seq2alnresidues <- 0
  for (i in 1:n) {
    chunkseq1aln <- substring(seq1aln, starts[i],
                              starts[i]+chunksize-1)
    chunkseq2aln <- substring(seq2aln, starts[i],
                              starts[i]+chunksize-1)
    # Find out how many gaps there are in chunkseq1aln:
    gaps1 <- countPattern("-",chunkseq1aln) # countPattern() is fromBiostrings library
    # Find out how many gaps there have been on this line of subject3:
    gaps2 <- countPattern("-",chunkseq2aln) # countPattern() is fromBiostrings library
    # Calculate how many residues of the first sequence we haveprinted so far in the alignment:
    seq1alnresidues <- seq1alnresidues + min(alnlen,chunksize) - gaps1
    # Calculate how many residues of the second sequence we haveprinted so far in the alignment:
    seq2alnresidues <- seq2alnresidues + min(alnlen,chunksize) - gaps2
    print(paste(chunkseq1aln,seq1alnresidues))
    print(paste(chunkseq2aln,seq2alnresidues))
    print(paste(' '))
  }
}


for (i in seq(1,length(a1$nam),2)) {
	
	b<-pairwiseAlignment(pattern=a1$seq[[i]],subject=a1$seq[[i+1]],type="global",gapOpening=10,gapExtension=0.5)
	s<-score(b)
	id<-pid(b)
	
        #r = BStringSet( c( toString( subject(b) ), toString( pattern(b)))) 
	l<-a1$nam[[i]]
	alignfile<-paste("./hpctmp/",l,"_print_align.txt",sep="")
	#writePairwiseAlignments(b,Matrix=NA,block.width=200)
	printPairwiseAlignment(b,100)
	write(id,file=paste("HPC_4N_B6XCAST_Alu_new_10k.txt.fa_200","_print_id.txt",sep=""),append="TRUE")

}

