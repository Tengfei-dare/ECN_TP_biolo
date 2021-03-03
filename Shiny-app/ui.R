#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#
library(shiny)
shinyUI(fluidPage(
  titlePanel("Tutorial Genetic Code"),
  tabsetPanel(
    # DNA sequence manipulation (length, GC-content, etc.)
    tabPanel("DNA sequence", 
             br(),
             h3("1 Load DNA sequence(s)"),
             fileInput("dna.file", "Fasta file (e.g. with the mysterious sequence)"),
             br(),
             h3("2 DNA content per sequence"),
             dataTableOutput("dna.content")
    ),
    # DNA sequence translation
    tabPanel("Protein sequence",
             br(),
             sidebarLayout(sidebarPanel(
               fileInput("dna.file.2", "Fasta file (e.g. with the mysterious sequence)"),
               br(),
               sliderInput("protein.length", "Minimum length for a protein", 1, 200, 80)
             ), mainPanel(
               fluidRow(column(width=12, wellPanel(
                 p("Compute proteins"), 
                 br(), 
                 uiOutput("aa.sequences")
               ))
               )
             ))
    ),
    # Protein sequence alignment
    tabPanel("Sequence alignment",
             br(),
             sidebarLayout(sidebarPanel(
               textAreaInput(inputId = "subject.seq",value = "MARVPRVRPPHGFALFLAKEEARKVKRLHGMLRSLLVYMLFLLVTLLASYGDASCHGHAYRLQSAIKQELH",
                             label = "Sequence to align against PKD1 of PKD2 references",
                             resize = 'both',
                             cols = 80,
                             rows = 10),
               selectInput("select.ref.seq", "Select reference sequence", c("PKD1_protein", "PKD2_protein"), "PKD1_protein", multiple = F),
               radioButtons("align.type", "Alignment type", c("global", "local"), "global"),
               selectInput("align.matrix", "Substitution matrix", c("BLOSUM45", "BLOSUM50",  "BLOSUM62", "BLOSUM80", "BLOSUM100", "PAM30", "PAM40", "PAM70", "PAM120", "PAM250"), "BLOSUM62"),
               sliderInput("gap.open.cost", "Gap opening cost", 0, 100, 10),
               sliderInput("gap.ext.cost", "Gap extension cost", 0, 100, 4)
             ),
             ## Exercise 7.4 - Fill here: sliderInput 2x
             mainPanel(h3("Alignment result"), 
                       verbatimTextOutput("align.out")
             ))
    )
  )
))