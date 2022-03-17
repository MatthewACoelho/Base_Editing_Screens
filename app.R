
#BE-viewer Finder Shiny app

#load libraries
library(shiny)
library(DT)
library(tidyverse)
library(ggrepel)
library(digest)
library(plotly)

data <- as_tibble(read.csv(file = "input.csv", header = TRUE)) %>%
    select(!(c(X, X.1))) %>%
    filter(is.na(zscore_proliferation) == FALSE)

pos <- position_jitter(seed = 2, width = 0.2, height = 0.2)

xy_str <- function(e) {
    paste0(round(e$x))
}

# Define UI for application
ui <- fluidPage(
    #logos
    img(height= 50, width = 160, src = "pic2.png"),
    img(height= 50, width = 160, src = "pic3.png"),
    br(),
    hr(),
    #title
    titlePanel("BE-view"),
    #description
    em("Interact with data from base editing mutagenesis screens to understand the functional consequence of variants in the IFNg signalling pathway.", 
         p("Proliferation (cell death) and FACS-based (MHC-I and PDL-1 expression) screening assays in HT-29 colorectal cancer cells.",
           p("CBE: cytidine base editor (C>T), ABE: adenine base editor (A>G)",
           )
         )
       ),
    br(),
    #image
    img(height= 80, width = 455, src = "pic.png"),
    br(),
    hr(),
    
    # Sidebar with input 
    sidebarLayout(
        sidebarPanel(
        selectInput("editor",
                        "editor", c("BE3-NGG (IFNG pathway)", "BE3-NGG (JAK1)", "BE4max-YE1-NGN (JAK1)", "BE3.9max-NGN (JAK1)", "ABE8e-NGN (JAK1)")),
        # Only show gene panel if the pathway screen selected
        conditionalPanel(
            condition = "input.editor == 'BE3-NGG (IFNG pathway)'",
            selectInput("gene",
                        "gene", c("JAK1", "STAT1", "IFNGR1", "IFNGR2", "SOCS1", "IRF1", "B2M", "JAK2")),
            ),
        # Only show screen panel if the available editor screens are selected
        conditionalPanel(
            condition = "input.editor == 'BE3-NGG (IFNG pathway)'| input.editor == 'BE3-NGG (JAK1)'",
            selectInput("screen",
                         "screen", c("proliferation", "FACS")),
        ),
        # Only show gene panel if the pathway screen selected
        conditionalPanel(
            condition = "input.editor == 'BE3-NGG (JAK1)' | input.editor == 'BE4max-YE1-NGN (JAK1)' | input.editor == 'BE3.9max-NGN (JAK1)' | input.editor == 'ABE8e-NGN (JAK1)'",
            selectInput("gene",
                        "gene", c("JAK1")),
        ),
        # Only show screen panel if the available editor screens are selected
        conditionalPanel(
            condition = "input.editor == 'BE4max-YE1-NGN (JAK1)' | input.editor == 'BE3.9max-NGN (JAK1)' | input.editor == 'ABE8e-NGN (JAK1)'",
            selectInput("screen",
                        "screen", c("proliferation")),
        ),
        selectInput("labels",
                    "labels", c("predicted amino acid change", "PTMs", "associated phenotypes", "citations")),
        hr(),
        downloadButton("downloadData", "download results"),
        width = 9,
        ),
        
        #MAIN
        mainPanel(
            #text
            hr(),
            textOutput("text"),
            
            #results
            plotlyOutput("plot"),
            br(),
            hr(),
            
            #data table
            DT::dataTableOutput("results"),
            
            #footnotes
            br(),
            hr(),
            h6("PTMs: post-translational modifications"),
            h6("gRNA off target summary: number of GRCh38 genomic positions with 0, 1, 2, 3, or 4 mismatches (please download data for this information)"),
            h6("guide: gRNA sequence (please download data for this information)"),
            
            #MIT license etc with footnotes
            hr(),
            #h5("please cite: X..."),
            a("publication", href = " "),
            br(),
            a("Matt Coelho's GitHub", href = "https://github.com/MatthewACoelho/"),
            br(),
            hr(),
            h6("License"),
            em("MIT: Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the \"Software\"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
            The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
            THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE."),
            hr(),
            h6("Shiny app development by Matt Coelho. Analysis by Matt Coelho. We used VEP for base editing predictions"),
            a("VEP", href = "https://www.ensembl.org/info/docs/tools/vep/index.html"),
            width = 9
        )
    )
)

# Define server logic required to make results table
server <- function(input, output) {
    
    text <- reactive({
        write <- paste0(input$gene, ", ", input$screen, ", ",input$editor)
    })
    
    results <- reactive({
                data <- data %>% 
                    filter(editor == input$editor) %>%
                    filter(Gene == input$gene) %>%
                    select(!c(guide, Amino_Acid_Position_simple, sgRNA_ID, off_target_summary_NGN, off_target_summary_NGG))
    })

    results_download <- reactive({
        data <- data %>% 
            filter(editor == input$editor) %>%
            filter(Gene == input$gene) %>%
            select(!c(Amino_Acid_Position_simple))
    })
    
   plot <- reactive({
    })
    
    
    #Render functions
    #text
    output$text <- renderText({
        text()
    })
    
    #plot
    output$plot <-renderPlotly({
        data <- data %>%
            filter(editor == input$editor) %>%
            filter(Gene == input$gene)
            
        if(input$screen == "proliferation" & input$labels == "predicted amino acid change")
            {plot_ly(data, 
                    x=~Amino_Acid_Position_simple, 
                    y=~zscore_proliferation,
                    color=~Consequence,
                    colors = c("black", "#fabd2f", "#cc241d", "blue"),
                    size = 2,
                    alpha = 0.8,
                    hoverinfo="text", 
                    text = ~paste(Amino_Acid_Change)
            ) %>% layout(xaxis=list(title = "amino acid position"), yaxis = list(title = "z-score"), legend = list(title = list(text = "consequence")))
        } 
        
        else if(input$screen == "proliferation" & input$labels == "associated phenotypes")
        {plot_ly(data, 
                 x=~Amino_Acid_Position_simple, 
                 y=~zscore_proliferation,
                 color=~Consequence,
                 colors = c("black", "#fabd2f", "#cc241d", "blue"),
                 size = 2,
                 alpha = 0.8,
                 hoverinfo="text", 
                 text = ~paste(Amino_Acid_Position, phenotype)
        ) %>% layout(xaxis=list(title = "amino acid position"), yaxis = list(title = "z-score"), legend = list(title = list(text = "consequence")))
        } 
        
        else if(input$screen == "proliferation" & input$labels == "PTMs")
        {plot_ly(data, 
                 x=~Amino_Acid_Position_simple, 
                 y=~zscore_proliferation,
                 color=~Consequence,
                 colors = c("black", "#fabd2f", "#cc241d", "blue"),
                 size = 2,
                 alpha = 0.8,
                 hoverinfo="text", 
                 text = ~paste(Amino_Acid_Position, PTM)
        ) %>% layout(xaxis=list(title = "amino acid position"), yaxis = list(title = "z-score"), legend = list(title = list(text = "consequence")))
        } 
        
        else if(input$screen == "proliferation" & input$labels == "citations")
        {plot_ly(data, 
                 x=~Amino_Acid_Position_simple, 
                 y=~zscore_proliferation,
                 color=~Consequence,
                 colors = c("black", "#fabd2f", "#cc241d", "blue"),
                 size = 2,
                 alpha = 0.8,
                 hoverinfo="text", 
                 text = ~paste(Amino_Acid_Position, Citation)
        ) %>% layout(xaxis=list(title = "amino acid position"), yaxis = list(title = "z-score"), legend = list(title = list(text = "consequence")))
        } 
        
        else if(input$screen == "FACS" & input$labels == "predicted amino acid change")
        {plot_ly(data, 
                 x=~Amino_Acid_Position_simple, 
                 y=~zscore_FACS,
                 color=~Consequence,
                 colors = c("black", "#fabd2f", "#cc241d", "blue"),
                 size = 2,
                 alpha = 0.8,
                 hoverinfo="text", 
                 text = ~paste(Amino_Acid_Change)
        ) %>% layout(xaxis=list(title = "amino acid position"), yaxis = list(title = "z-score"), legend = list(title = list(text = "consequence")))
        } 
        
        else if(input$screen == "FACS" & input$labels == "associated phenotypes")
        {plot_ly(data, 
                 x=~Amino_Acid_Position_simple, 
                 y=~zscore_FACS,
                 color=~Consequence,
                 colors = c("black", "#fabd2f", "#cc241d", "blue"),
                 size = 2,
                 alpha = 0.8,
                 hoverinfo="text", 
                 text = ~paste(Amino_Acid_Position, phenotype)
        ) %>% layout(xaxis=list(title = "amino acid position"), yaxis = list(title = "z-score"), legend = list(title = list(text = "consequence")))
        } 
        
        else if(input$screen == "FACS" & input$labels == "PTMs")
        {plot_ly(data, 
                 x=~Amino_Acid_Position_simple, 
                 y=~zscore_FACS,
                 color=~Consequence,
                 colors = c("black", "#fabd2f", "#cc241d", "blue"),
                 size = 2,
                 alpha = 0.8,
                 hoverinfo="text", 
                 text = ~paste(Amino_Acid_Position, PTM)
        ) %>% layout(xaxis=list(title = "amino acid position"), yaxis = list(title = "z-score"), legend = list(title = list(text = "consequence")))
        } 
        
        else if(input$screen == "FACS" & input$labels == "citations")
        {plot_ly(data, 
                 x=~Amino_Acid_Position_simple, 
                 y=~zscore_FACS,
                 color=~Consequence,
                 colors = c("black", "#fabd2f", "#cc241d", "blue"),
                 size = 2,
                 alpha = 0.8,
                 hoverinfo="text", 
                 text = ~paste(Amino_Acid_Position, Citation)
        ) %>% layout(xaxis=list(title = "amino acid position"), yaxis = list(title = "z-score"), legend = list(title = list(text = "consequence")))
        } 
    })
    
    #hover
    output$info <- renderText({
        info()
    })
    
    #show results data table in the shiny app
    output$results <- DT::renderDataTable({
        results()
    })
    
    # Downloadable txt of data
    output$downloadData <- downloadHandler(
        filename = function() {
            paste(c(input$gene, "_mutagenesis.txt"), collapse = "")
        },
        content = function(file) {
            write.table(results_download(), file, sep ="\t")
        }
    )
}

# Run the application 
shinyApp(ui = ui, server = server)
