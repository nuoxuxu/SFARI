library(bslib)

ui <- page_fillable(
  selectizeInput( 
    "select_genes", 
    "Select genes from list below:", 
    choices = c("SCN2A")
  ),
  plotOutput("ggtranscript_plot"),
)