library(bslib)

ui <- page_fluid(
  selectizeInput( 
    "select_genes", 
    "Select genes from list below:", 
    choices = NULL
  ),
  plotOutput("ggtranscript_plot"),
)