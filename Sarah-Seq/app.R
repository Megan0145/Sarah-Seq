#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(vcfR)
library(tidyr)
library(ggplot2)
library(scales)
library(VariantAnnotation)
library(matrixStats)
library(karyoploteR)
# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Sarah-Seq"),

    # Sidebar with a slider input for number of bins 
    tabsetPanel(
        tabPanel("Upload File",
            sidebarLayout(
                sidebarPanel(
                    h4("Upload a File"),
                    fileInput("file1", "Choose VCF file",
                      multiple = TRUE,
                      accept = c(".vcf")),
                    tags$hr(),
                    radioButtons("assem", "Choose an assembly", c('hg38', 'hg19'), selected = 'hg38'),
                    tags$hr()
        )
        , 
                mainPanel(
                    verbatimTextOutput("success", placeholder = F))
    )),
        tabPanel("QC",
             sidebarLayout(
                 sidebarPanel(
                     selectInput("chr", "Choose a chromosome", c('All', seq(1:22), 'X', 'Y'), selected = 'All' ),
                     tags$hr(),
                     h5(tags$b("Number of Variants")),
                     verbatimTextOutput("n_vars", placeholder = TRUE)
                 ),
                 mainPanel(
                    plotOutput("QCPlot")
                   ))),
        tabPanel("Visualise Distribution of Variants",
            #sidebarLayout(
                #sidebarPanel(
                    #selectInput("chr", "Choose a chromosome", c('All', seq(1:22), 'X', 'Y'), selected = 'All' ),
                    #tags$hr(),
                    #h5(tags$b("Number of Variants")),
                    #verbatimTextOutput("n_vars", placeholder = TRUE)
                    #),
                #mainPanel(
                    plotOutput("karyoPlot")
            #))
            )
))

# Define server logic required to draw a histogram
server <- function(input, output) {
    options(shiny.maxRequestSize=500*1024^2)
    vcfInput <- reactive({
        req(input$file1)
        vcf    <- readVcf(input$file1$datapath)
        
    })
    dataInput <- reactive({
        if(input$chr == 'All'){
            data <- vcfInput()
        } 
        else{
            data <- subset(vcfInput(), vcfInput()@rowRanges@seqnames == paste('chr', input$chr, sep = ""))
        }
    })
   
    output$success <- renderText({
        req(input$file1)
        print(paste("File with ", nrow(vcfInput()), " variants from ", rownames(vcfInput()@colData), " successfully uploaded!"))
    })    
    
    output$karyoPlot <- renderPlot({
        if (input$chr == 'All'){
            kp <- plotKaryotype(input$assem)
            kpPlotDensity(kp, data=rowRanges(dataInput()), col = '#9fe0d3')
            }
        else{
            kp <- plotKaryotype(input$assem, chromosomes = paste('chr',input$chr, sep=""))
            kp <- kpPlotDensity(kp, data=rowRanges(dataInput()), col = 'indianred')
            kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, cex=0.8)
            kpAbline(kp, h=mean(kp$latest.plot$computed.values$density), lty=2, ymax=kp$latest.plot$computed.values$max.density, r0=0.68, r1=1)
            kpAddBaseNumbers(kp)
            title(main = paste('Distribution of variants on chromosome', input$chr))
            }
    })
    output$QCPlot <- renderPlot({
        qc_info <- data.frame(QUAL = dataInput()@fixed$QUAL, DP = dataInput()@info$DP, 
                              MQ = dataInput()@info$MQ, stringsAsFactors = F)
        qc_info <- gather(qc_info, measure)
        ggplot(qc_info, aes(x = value, col = measure, fill = measure)) +
            geom_histogram(show.legend = F, alpha = 0.4) +
            facet_wrap(~measure, scales = 'free') +
            labs(title = 'Distribution of QC metrics', x = 'Score', y = 'Count') +
            theme_linedraw() +
            theme(axis.text.x = element_text(angle = 90))
    })
    
    output$n_vars <- renderText({
            nrow(dataInput())
        })
}

# Run the application 
shinyApp(ui = ui, server = server)
