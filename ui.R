library(shiny)

scheme.list <- c("20/20" = "2020", "NNK/S" = "NNK", "NNB" = "NNB", "NNN" = "NNN", "20/20 (-C)" = "2020C", "NNK/S (-C)" = "NNKC", "NNB (-C)" = "NNBC", "NNN (-C)" = "NNNC")

makePopover <- function(popover, direction = "right") {
    return(HTML(paste("<i class='fa fa-info-circle' data-toggle='popover' data-placement='", direction, "' rel = 'popover' data-content='", paste0(popover, "<br><br>___________________________________________________________"), "'></i>", sep = "")))
}

source("R/strings.R")

shinyUI(fluidPage(
    
    titlePanel("", windowTitle = "PeLiCa - Peptide Library Calculator"),
    
    sidebarLayout(
        sidebarPanel(
            img(src = "img/pelica.gif", width = 350),
            
            hr(),
            
            includeHTML("www/html/additional.html"),
            includeCSS("www/css/additional.css"),
            includeScript("www/js/additional.js"),
            
            tags$head(includeScript("www/js/google-analytics.js"), tags$link(rel="shortcut icon", href="img/favicon.ico")),
            
            navbarPage("Configuration",
                       tabPanel("Basic",
                                
                                h4("Library"),
                                
                                fluidRow(
                                    column(width = 7, 
                                           selectizeInput("scheme", label = "Scheme", choices = c(scheme.list, "Custom"))
                                    ), 
                                    column(width = 1, offset = 3, 
                                           makePopover(tooltip_scheme)
                                    )
                                ),
                                
                                conditionalPanel("input.scheme == 'Custom'",
                                                 fluidRow(
                                                     column(width = 9, 
                                                            fluidRow(
                                                                column(width = 8,
                                                                       textInput("custom_aacid1", "Amino Acids", value="ACDEFGHIKLMNPQRSTVWY"),
                                                                       uiOutput("custom_aacid"),
                                                                       textInput("custom_aacid0", "", value="*"),
                                                                       actionButton("custom_new", "Add New AA Class")
                                                                ),
                                                                
                                                                column(width = 4,
                                                                       numericInput("custom_c1", "Codons", value = 1),
                                                                       uiOutput("custom_c"),
                                                                       numericInput("custom_c0", "", value=0),
                                                                       actionButton("custom_done", "Done")
                                                                )
                                                            )
                                                     ), 
                                                     column(width = 1, offset = 1, 
                                                            makePopover(tooltip_custom)
                                                     )
                                                 ),
                                                 
                                                 hr()
                                ),
                                
                                fluidRow(
                                    column(width = 4, 
                                           numericInput("size_x", "Size (X)", min = 1.0, max = 9.9, step = 0.1, value = 1)
                                    ),
                                    column(width = 3,
                                           helpText("x 10^") 
                                    ),
                                    column(width = 3,
                                           numericInput("size_y", "Size (Y)", min = 1, max = 25, step = 1, value = 8)
                                    ),
                                    column(width = 1, 
                                           makePopover(tooltip_size)
                                    )
                                ),
                                
                                hr(),
                                
                                h4("Peptide"),
                                
                                fluidRow(
                                    column(width = 7, 
                                           sliderInput("length", "Length", min = 1, max = 20, step = 1, value = 7)
                                    ),
                                    column(width = 1, offset = 3, 
                                           makePopover(tooltip_length)
                                    )
                                ),
                                
                                fluidRow(
                                    column(width = 7, 
                                           textInput("peptide", "Sequence", value = "HENNING")
                                    ),
                                    column(width = 1, offset = 3, 
                                           makePopover(tooltip_peptide)
                                    )
                                )
                       ),
                       
                       tabPanel("Advanced",
                                h4("Other"),
                                
                                fluidRow(
                                    column(width = 7, 
                                           numericInput("digits", "Digits", value = 4, min = 0, max = 8)
                                    ),
                                    column(width = 1, offset = 3, 
                                           makePopover(tooltip_digits)
                                    )
                                ),
                                
                                fluidRow(
                                    column(width = 7, 
                                           selectizeInput("schemes", label = "Display Schemes", choices = scheme.list, multiple = TRUE)
                                    ),
                                    column(width = 1, offset = 3, 
                                           makePopover(tooltip_schemes)
                                    )
                                ),
                                
                                fluidRow(
                                    column(width = 7, 
                                           selectizeInput("lengths", label = "Display Lengths", choices = 1:20, multiple = TRUE)
                                    ),
                                    column(width = 1, offset = 3, 
                                           makePopover(tooltip_lengths)
                                    )
                                )
                       )
            )
        ),
        
        mainPanel(
            navbarPage("Results",
                       tabPanel("Welcome",
                                fluidRow(
                                    column(width = 11, 
                                           h3("PeLiCa - The Peptide Library Calculator")
                                    ),
                                    column(width = 1, 
                                           makePopover(tooltip_welcome, direction = "bottom")
                                    )
                                ),
                                
                                hr(),
                                
                                h4("About"),
                                HTML(text_welcome),
                                
                                hr(),
                                
                                h4("Authors"),
                                HTML(text_authors)
                       ),
                       
                       tabPanel("Summary",
                                fluidRow(
                                    column(width = 11, 
                                           h4("Library Properties")                             ),
                                    column(width = 1, 
                                           makePopover(tooltip_summary_library, direction = "bottom")
                                    )
                                ),
                                dataTableOutput("properties_table"),
                                
                                hr(),
                                
                                fluidRow(
                                    column(width = 11, 
                                           h4("Peptide Properties")                             ),
                                    column(width = 1, 
                                           makePopover(tooltip_summary_peptide, direction = "bottom")
                                    )
                                ),
                                dataTableOutput("pep_properties_table"),
                                
                                hr(),
                                
                                fluidRow(
                                    column(width = 11, 
                                           h4("Scheme Properties")                             ),
                                    column(width = 1, 
                                           makePopover(tooltip_summary_scheme, direction = "bottom")
                                    )
                                ),
                                dataTableOutput("scheme_table"),
                                
                                hr(),
                                
                                fluidRow(
                                    column(width = 11, 
                                           h4("Peptide Sample")                             ),
                                    column(width = 1, 
                                           makePopover(tooltip_summary_sample, direction = "bottom")
                                    )
                                ),
                                dataTableOutput("sample_table")
                       ),
                       
                       tabPanel("Inclusion",
                                fluidRow(
                                    column(width = 11, 
                                           tags$b(textOutput("inclusion_text"))
                                    ),
                                    column(width = 1, 
                                           makePopover(tooltip_inclusion, direction = "bottom")
                                    )
                                ),
                                fluidRow(
                                    column(width = 12,                                      
                                           tags$b(textOutput("inclusion_text_pep")),
                                           
                                           hr(),
                                           
                                           textOutput("inclusion_caption"),
                                           plotOutput("inclusion_plot"),
                                           
                                           hr(),
                                           
                                           textOutput("inclusion_caption_schemes"),
                                           plotOutput("inclusion_plot_schemes")
                                    )
                                )
                       ),  
                       tabPanel("Neighborhood",
                                fluidRow(
                                    column(width = 11, 
                                           tags$b(textOutput("neighborhood_text_degree1")),
                                           tags$b(textOutput("neighborhood_text_degree2")),       
                                           
                                           br(),
                                           
                                           tags$b(textOutput("neighborhood_text_mine_degree1")),
                                           tags$b(textOutput("neighborhood_text_mine_degree2"))
                                    ),
                                    column(width = 1, 
                                           makePopover(tooltip_neighborhood, direction = "bottom")
                                    )
                                ),
                                
                                hr(),
                                
                                fluidRow(
                                    column(width = 12,          
                                           textOutput("neighborhood_caption"),
                                           plotOutput("neighborhood_plot"),
                                           
                                           hr(),
                                           
                                           textOutput("neighborhood_caption_mine"),
                                           plotOutput("neighborhood_plot_mine")
                                    )
                                )
                       ),
                       
                       tabPanel("Coverage",
                                fluidRow(
                                    column(width = 11, 
                                           tags$b(textOutput("coverage_text"))
                                    ),
                                    column(width = 1, 
                                           makePopover(tooltip_coverage, direction = "bottom")
                                    )
                                ),
                                
                                hr(),
                                
                                fluidRow(
                                    column(width = 12,                                        
                                           textOutput("coverage_caption"),
                                           plotOutput("coverage_plot")
                                    )
                                )
                       ),
                       
                       tabPanel("Efficiency",
                                fluidRow(
                                    column(width = 11, 
                                           tags$b(textOutput("efficiency_text"))
                                    ),
                                    column(width = 1, 
                                           makePopover(tooltip_efficiency, direction = "bottom")
                                    )
                                ),
                                
                                hr(),
                                
                                fluidRow(
                                    column(width = 12,                                        
                                           textOutput("efficiency_caption"),
                                           plotOutput("efficiency_plot")
                                    )
                                )
                       ),
                       
                       tabPanel("Diversity",
                                fluidRow(
                                    column(width = 11,         
                                           tags$b(textOutput("pepdiversity_text")),
                                           textOutput("pepdiversity_caption")
                                    ),
                                    column(width = 1, 
                                           makePopover(tooltip_diversity, direction = "bottom")
                                    )
                                ),
                                
                                hr(),
                                
                                fluidRow(
                                    column(width = 12, 
                                           tags$b(textOutput("diversity_text")),
                                           textOutput("diversity_caption"),
                                           tableOutput("diversity_table")
                                    )
                                )
                       )
            )
        )
    )
))
