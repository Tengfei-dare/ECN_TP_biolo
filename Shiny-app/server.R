#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#
library(shiny)
library(seqinr)
library(Biostrings)

# 'util' file for genetic code logic and strings cleaning
source("util.R")

# server logic to read fasta files
shinyServer(function(input, output) {
  
  output$dna.content <- renderDataTable({
    
    # save inFile for internal logic
    inFile <- input$dna.file
    
    if (is.null(inFile))
      return(NULL)
    
    dna.sequences <- read.fasta(inFile$datapath, seqtype = "DNA", forceDNAtolower = T)
    
    dna.content <- compute.dna.metrics(dna.sequences)
    
    return(dna.content)
  })
  
  output$aa.sequences <- renderUI({
    inFile.2 <- input$dna.file.2
    p.protein <- input$protein.length #get the slide input
    
    if (is.null(inFile.2))
      return(NULL)
    # load fasta file - Exercise 6.4
    dna.sequences <- read.fasta(inFile.2$datapath, seqtype = "DNA", forceDNAtolower = T)
      # find mysterious proteins - Exercise 6.4
      list.of.proteins <- find.mysterious.proteins(dna.sequences$seq0, p.protein)
      
      HTML(proteins.to.html(list.of.proteins))
})
  
  output$align.out <- renderPrint({
    
    if (input$select.ref.seq == "PKD2_protein"){
      ref.sequence <- readAAStringSet("../TP-Student/Students/Reference_protein_seq.txt")[2]
    }
    else{
      ref.sequence <- readAAStringSet("../TP-Student/Students/Reference_protein_seq.txt")[1]
    }
    
      in.seq <- clean.sequence(input$subject.seq)
    if(is.null(in.seq))
      return(NULL)
    pair.alignment <- pairwiseAlignment(pattern = ref.sequence,
                                        subject = in.seq,
                                        substitutionMatrix = input$align.matrix,
                                        type = input$align.type,
                                        gapOpening = input$gap.open.cost,
                                        gapExtension = input$gap.ext.cost)
      
      return(pair.alignment)
  })
  
})
