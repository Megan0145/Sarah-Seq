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
    h1(id ='title', "Sarah-Seq"),
    tags$style(HTML("@import url('//fonts.googleapis.com/css?family=Lobster|Cabin:400,700');
                    #title{color: #3588b8;
                    font-family: 'Lobster', cursive;}")),

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
                  fluidRow(align = 'center',
                           style = "border: 8px double white;",
                           img(src='michael.gif', align = "center", border = "5cm"),
                           tags$hr(),
                           verbatimTextOutput("success", placeholder = F)))
    ))
    ,
        tabPanel("QC",
             sidebarLayout(
                 sidebarPanel(
                     selectInput("chr", "Choose a chromosome", c('All', seq(1:22), 'X', 'Y'), selected = 'All' ),
                     h5(tags$b("Number of variants")),
                     verbatimTextOutput("n_vars", placeholder = TRUE),
                     tags$hr(),
                     h4(tags$b("Filtering")),
                     h5("Include variants that have:"),
                     sliderInput("mq_fil", "Mapping quality scores between", min = 1, max = 60, value = c(1,60), step = 1),
                     sliderInput("dp_fil", "Read depth scores between", min = 1, max = 7000, value = c(1,7000), step = 1),
                     sliderInput("qual_fil", "Quality scores between", min = 1, max = 250000, value = c(1, 250000)),
                     tags$hr(),
                     fluidRow(align = 'center',
                     actionButton("filter", tags$b("Update View")))
                 ),
                 mainPanel(
                   h4(tags$b("QC statistics:")),
                   tableOutput("qc_stats"),
                   tags$hr(),
                   h4(tags$b("QC plots:")),
                   plotOutput("QCPlot")
                   ))),
        tabPanel("Visualise Distribution of Variants",
                    plotOutput("karyoPlot")),
        tabPanel("PRS"),
        tabPanel("Ancestry"),
        tabPanel("Mitochondrial Genome")
            
))

# Define server logic required to draw a histogram
server <- function(input, output) {
    options(shiny.maxRequestSize=500*1024^2)
    vcfInput <- reactive({
        req(input$file1)
        vcf    <- readVcf(input$file1$datapath)
        
    })
    dataInput <- eventReactive(
        input$filter,{
        if(input$chr == 'All'){
            predata <- vcfInput()
        } 
        else{
            predata <- subset(vcfInput(), vcfInput()@rowRanges@seqnames == paste('chr', input$chr, sep = ""))
            
        }
      data <- subset(predata, predata@info$MQ >= input$mq_fil[1] & predata@info$MQ <= input$mq_fil[2]
                              & predata@info$DP >= input$dp_fil[1] & predata@info$DP <= input$dp_fil[2]
                              & predata@fixed$QUAL >= input$qual_fil[1] & predata@fixed$QUAL <= input$qual_fil[2])
    }, ignoreNULL = F)
    
    output$success <- renderText({
        req(input$file1)
        print(paste("File with ",nrow(vcfInput()), " variants from ", rownames(vcfInput()@colData), " successfully uploaded!"))
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
    qcInfo <- reactive({qc_info <- data.frame(DP = dataInput()@info$DP, MQ = dataInput()@info$MQ,
                                              QUAL = dataInput()@fixed$QUAL,  stringsAsFactors = F)})
    output$QCPlot <- renderPlot({
        qc_info <- gather(qcInfo(), measure)
        ggplot(qc_info, aes(x = value, col = measure, fill = measure)) +
            geom_histogram(show.legend = F, alpha = 0.4, bins = 50) +
            facet_wrap(~measure, scales = 'free') +
            labs(title = 'Distribution of QC metrics', x = 'Score', y = 'Count') +
            theme_linedraw() +
            theme(axis.text.x = element_text(angle = 90))
    })
    
    output$n_vars <- renderText({
            nrow(dataInput())
        })
    output$qc_stats <- renderTable(rownames = T,{
      do.call(cbind, lapply(qcInfo(), summary)) })
}

# Run the application 
shinyApp(ui = ui, server = server)
