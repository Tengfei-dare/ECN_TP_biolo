# TP R - Code g¨¦n¨¦tique
#
# Aubin Dauny-Tengfei HAN
# groupe 1-8
#

############
# packages #
############
library(seqinr)
library(shiny)
library(micropan)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")
library(Biostrings)

##############
# load files #
##############
sequences <- read.fasta(file = "./TP-Student/Students/Mysterious_seq.txt", 
                        seqtype = "DNA",
                        forceDNAtolower = T)



################
# Exercice 2.1 #
# returns the number of occurrences
compte <- function(letter, seq){
  occur = seq == letter
  letter.occur <- c(letter, sum(occur))
  # to print letter and the number for a better view
  print(letter.occur)
  # return(sum(occur))
}

sequence <- c("a","a","t","g","a","g","c","t","a","g","c","t","g") # load seq
compte("a", sequence) # compte how many a

# compte each letter in seq
for (elt in unique(sequences$seq0)){
  compte(elt, sequences$seq0)
}

################
# Exercice 3.1 #
# options
# ?read.fasta

# Print the first 50 nucleotides
head(sequences$seq0, n = 50)
# Print the length
length(sequences$seq0)
# Print the number of occurences
table(sequences$seq0)

################
# Exercice 3.2 #
# print the proportion of each nucl¨¦otide
proportions <- table(sequences$seq0)/length(sequences$seq0)
# GC-content = 
proportions[[2]]+proportions[[3]]

################
# Exercice 3.3 #
complementary <- function(seq){
  com.seq <- chartr("A", "t", chartr("t", "a", 
                          chartr("a", "A", chartr("C", "g", 
                                                  chartr("g", "c", chartr("c", "C", seq))))))
  return(com.seq)
}
# print the complementary sequence
complementary(sequences$seq0)
# print the reverse complementary
rev(complementary(sequences$seq0))


################
# Exercice 5.1 #
transcribe <- function(cDNA){
  RNA <- chartr("A", "u", chartr("t", "a", 
                                 chartr("a", "A", 
                                        chartr("C", "g", 
                                               chartr("g", "c", 
                                                      chartr("c", "C", 
                                                             cDNA))))))
  return(RNA)
}
transcribe(sequences$seq0)


################
# Exercice 6.2 #
# read codon table and add column names
codon.table <- read.table(file = "./TP-Student/Code_shortcut/Genetic_code.txt", 
                          col.names = c("codon", "aa", "letter"),
                          stringsAsFactors = F)
codon.table$codon <- tolower(codon.table$codon) # use only lower case characters
rownames(codon.table) <- codon.table$codon # use codons as table keys
head(codon.table)

# translates a DNA sequence into an amino acid chain
dna2peptide <- function(sequence, frame=0, sens='F'){
  # adapted from seqinr
  if(sens=='R') {
    sequence <- rev(complementary(sequence))
  }
  # length of the DNA sequence, %/% stands for euclidean division
  l <- 3 * ((length(sequence) - frame)%/%3)
  # list of indexes for first nucleotides in codons
  c1 <- seq(frame+1, frame+l, 3) # from, to, by
  # list of codons
  codons <- paste0(sequence[c1], sequence[c1+1], sequence[c1+2])
  # translate codons to amino acids
  peptide <- codon.table[codons, "letter"]
  return(peptide)
}
# test the function
# dna2peptide(c('a', 't', 'g', 't', 't', 'c', 't', 't', 't', 'a'), 0, 'F')
dna2peptide(sequences$seq0)

################
# Exercice 6.3 #
# takes a cDNA sequence as input and returns a list of proteins
find.mysterious.proteins <- function(sequence){
  out = list()
  for (f in 0:2){ # vary the frame
    for (s in c("F", "R")){ # vary the direction
      # gregexpr: find all the matches with a list
      temp = gregexpr('M[ACDEFGHIKLMNPQRSTVWY]{80,}(X|$)',
                      paste0(dna2peptide(sequence, f, s), collapse = ""),
                      extract = TRUE)
      if (temp != ""){ # if find, add to list
        for (elt in temp){
          out = append(out, elt)
        }
      }
    }
  }
  return(out)
}

mysterious.proteins <- find.mysterious.proteins(sequences$seq0)
head(mysterious.proteins)

#save the result in a file
write.table(mysterious.proteins, 
            file = "./TP-Student/Students/mysterious_proteins.txt",
            quote = FALSE)


################
# Exercice 7.1 #
ref.prot.sequences <- readAAStringSet('./TP-Student/Students/Reference_protein_seq.txt')
target.prot.sequences <- readAAStringSet('./TP-Student/Students/Protein_sequences.txt')
summary(ref.prot.sequences)
summary(target.prot.sequences)

# align mysterious proteins
my.sequence <- mysterious.proteins[[1]] # the longest translated amino acid sequence
pair.alignment <- pairwiseAlignment(pattern = ref.prot.sequences,
                                    subject = my.sequence,
                                    substitutionMatrix = "BLOSUM62",
                                    type = "global")
writePairwiseAlignments(pair.alignment, file="./TP-Student/Students/alignment.my.sequence.PKD1-2.txt")

# align target proteins
target.pair.alignment <- pairwiseAlignment(pattern = ref.prot.sequences,
                                           subject = target.prot.sequences$seq8,
                                           substitutionMatrix = "BLOSUM62",
                                           type = "global")
writePairwiseAlignments(target.pair.alignment, file="./TP-Student/Students/alignment.target.sequence.PKD1-2.txt")


################
# Exercice 7.2 #
protein32seq <- c()
for (i in 1:32) {
  protein1seq <- pairwiseAlignment(pattern = target.prot.sequences,
                                   subject = target.prot.sequences[[i]],
                                   substitutionMatrix = "BLOSUM62",
                                   type = "global")@score
  protein32seq <- rbind(protein32seq, protein1seq)
}
scores <- as.matrix(protein32seq)
row.names(scores) <- 1:32

heatmap(scores, Rowv = NA, Colv = NA, symm = T, 
        main = "Alignment scores heatmap for 32 protein sequences")

# Interchanging columns and lines to see the two clusters distinctly
map <- scores
map2 = cbind(map[,1:2], map[,4:7], map[,11:17], map[,20:21], map[,24], map[,26], map[,29], map[,3], map[,8:10], map[,18:19], map[,22:23], map[,25], map[,27:28], map[,30:32])
map3 = rbind(map2[1:2,], map2[4:7,], map2[11:17,], map2[20:21,], map2[24,], map2[26,], map2[29,], map2[3,], map2[8:10,], map2[18:19,], map2[22:23,], map2[25,], map2[27:28,], map2[30:32,])
map4 = rbind(map3[1:6,], map3[8:12,], map3[15:17,], map3[7,], map3[13,], map3[14,], map3[18,], map3[19:32,])
map5 = cbind(map4[,1:6], map4[,8:12], map4[,15:17], map4[,7], map4[,13], map4[,14], map4[,18], map4[,19:32])
heatmap(map5, Rowv = NA, Colv = NA, symm = T, main = "Alignment scores heatmap for 32 protein sequences")
map5
