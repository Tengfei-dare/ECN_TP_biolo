compute.dna.metrics <- function(dna.sequences){
  # Bonus exercise 7.5: combine the 6 sapply calls
  dna.content <- data.frame(names=names(dna.sequences))
  
  dna.content$a <- sapply(dna.sequences, compte, 'a')
  dna.content$c <- sapply(dna.sequences, compte, 'c')
  dna.content$g <- sapply(dna.sequences, compte, 'g')
  dna.content$t <- sapply(dna.sequences, compte, 't')
  
  dna.content$length <- sapply(dna.sequences, length)
  dna.content$GC.content <- sapply(dna.sequences, GC)
  
  return(dna.content)
}

compte <- function(seq, letter){
  # Exercise 4.1 with the Exercise 2.1 function
  occur = seq == letter
  letter.occur <- sum(occur)
  return(letter.occur)
}

codon.table <- read.table(file = "../TP-Student/Code_shortcut/Genetic_code.txt", 
                          col.names = c("codon", "aa", "letter"),
                          stringsAsFactors = F)
codon.table$codon <- tolower(codon.table$codon) # use only lower case characters
rownames(codon.table) <- codon.table$codon # use codons as table keys

complementary <- function(s){
  return(chartr("A", "t", chartr("t", "a", chartr("a", "A", chartr("C", "g", chartr("g", "c", chartr("c", "C", s)))))))
}

dna2peptide <- function(sequence, frame=0, sens='F'){
  # Exercise 6.4
  if(sens=='R') {
    sequence <- rev(complementary(sequences$seq0))
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

find.mysterious.proteins = function(sequence, p.length){
  out = list()
  for (f in 0:2){
    for (s in c("F", "R")){
      temp = gregexpr(paste0('M[ACDEFGHIKLMNPQRSTVWY]{',p.length ,',}(X|$)', collapse = ""), paste0(dna2peptide(sequence, f, s), collapse = ""), extract = TRUE)
      if (temp != ""){
        for (elt in temp){
          out = append(out, elt)
        }
      }
    }
  }
  return(out)
}

clean.sequence <- function(dirty.seq){
  return(gsub('[^A-Z]', replacement = "", toupper(dirty.seq)))
}

proteins.to.html <- function(list.of.proteins){
  # Every 60 amino acids, insert a new line (html code)
  proteins.cut <- gsub("(.{61,}?)", "\\1</br>", list.of.proteins)
  # Use console typo in html
  proteins.in.span <- paste0("<p style='font-family:monospace'>",proteins.cut,"</span>")
  return(paste0(proteins.in.span, collapse = "</br></br>"))
}

# more to come