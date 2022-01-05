
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
    hr(),
    #title
    titlePanel("BE-viewer"),
    hr(),
    #description
    em("interact with base editing mutagenesis data - response to IFNg in HT-29 colorectal cancer cells", 
         p("we used proliferation (cell death) and FACS-based (MHC-I and PDL-1 expression) screening assays",
           p("CBE: cytosine base editor, ABE: adenine base editor"
           )
         )
       ),
    br(),
    
    # Sidebar with input 
    sidebarLayout(
        sidebarPanel(
        selectInput("gene",
                  "gene", c("JAK1", "STAT1", "IFNGR1", "IFNGR2", "SOCS1", "IRF1", "B2M", "JAK2")),
        selectInput("screen",
                  "screen", c("proliferation", "FACS")),
        selectInput("editor",
                    "editor", c("BE3_NGG_JAK1_focused", "BE3_NGG_IFNG_pathway", "BE4max_YE1_NGN", "BE3.9max_NGN", "ABE8e_NGN")),
        hr(),
        downloadButton("downloadData", "download results"),
        ),
        
        #MAIN
        mainPanel(
            br(),
            #image
            img(height= 85, width = 455, src = "pic.png"),
            br(),
            hr(),
            
            #text
            textOutput("text"),
            
            #results
            plotlyOutput("plot"),
            br(),
            hr(),
            
            #data table
            DT::dataTableOutput("results"),
            
            #MIT license etc with footnotes
            br(),
            hr(),
            #h5("please cite: X..."),
            a("publication", href = " "),
            br(),
            a("Matt Coelho's GitHub", href = "https://github.com/MatthewACoelho/"),
            br(),
            hr(),
            h5("License"),
            h6("Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the \"Software\"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
            The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
            THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE."),
            hr(),
            h6("Shiny app development by Matt Coelho. Analysis by Matt Coelho. We used VEP for base editing predictions"),
            a("VEP", href = "https://www.ensembl.org/info/docs/tools/vep/index.html"),
         
        )
    )
)

# Define server logic required to make results table
server <- function(input, output) {
    
    text <- reactive({
        write <- paste0(input$gene, ", ", input$screen, ", ",input$editor)
    })
    
    results <- reactive({
                data <- data %>% filter(Gene == input$gene) %>%
                    filter(editor == input$editor) %>%
                    select(!c(Amino_Acid_Position_simple))
    })
    
    
    
    #Render functions

    #text
    output$text <- renderText({
        text()
    })
    
    #plot
    output$plot <- renderPlotly({
        data <- data %>% filter(Gene == input$gene) %>%
            filter(editor == input$editor)
        
        if (input$screen == "proliferation") {
            plot_ly(data, 
                    x=~Amino_Acid_Position_simple, 
                    y=~zscore_proliferation,
                    color=~Consequence,
                    colors = c("black", "#fabd2f", "#cc241d", "blue"),
                    size = 2,
                    alpha = 0.8,
                    hoverinfo="text", 
                    text = ~paste0(Amino_Acid_Change)
            ) %>% layout(xaxis=list(title = "amino acid position"), yaxis = list(title = "z-score"), legend = list(title = list(text = "consequence")))
        } else {
            plot_ly(data, 
                    x=~Amino_Acid_Position_simple, 
                    y=~zscore_FACS,
                    color=~Consequence,
                    colors = c("black", "#fabd2f", "#cc241d", "blue"),
                    size = 2,
                    alpha = 0.8,
                    hoverinfo="text",
                    text = ~paste0(Amino_Acid_Change)
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
            write.table(results(), file, sep ="\t")
        }
    )
}

# Run the application 
shinyApp(ui = ui, server = server)
