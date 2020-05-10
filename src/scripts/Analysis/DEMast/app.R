#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("LncRNAs DE overview"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("fdr",
                        "FDR Threshold:",
                        min = 0,
                        max = 1,
                        value = 0.1),
            sliderInput("fc",
                            "Fold Change Threshold:",
                            min = 1,
                            max = 3,
                            value = 1, step = 0.1)
            
        ),
    


        # Show a plot of the generated distribution
        mainPanel(
          
           plotOutput("infectedvsbystanders"),
           #plotOutput("infectedvsbystanders_lnc"), 
           plotOutput("celltypes")
        )
    )
)
library(plotly)
#hurdle_summary_celltypes_lnc <- readRDS("hurdle_summary_celltypes_lnc.rds")
# Define server logic required to draw a histogram
server <- function(input, output) {

    output$infectedvsbystanders <- renderPlot({
        # generate bins based on input$bins from ui.R
        # draw the histogram with the specified number of bins
        FDRTHRESHOLD <- input$fdr
        FCTHRESHOLD <- log2(input$fc)
        df_filtered <- fcHurdle_all_infectedVSbystanders[fdr<FDRTHRESHOLD & abs(coef)>FCTHRESHOLD,]
        df_plot_ready <- df_filtered %>% group_by(type) %>% dplyr::mutate(count= n_distinct(primerid))
        df_plot_ready[df_plot_ready$type == "lncRNA",]$primerid
        ggplot(df_plot_ready, aes(x = type, y =count, fill = type))+
            geom_bar(stat = "identity", color = "black", position =position_dodge())+
            theme_classic()+
            scale_fill_brewer(palette = "Set2")+
            ylim(0,as.integer(max(df_plot_ready$count)))+
            labs(x = "", y = "# DE Genes ", title = "Infected vs Bystanders in Myeloid 24H")
    })
    
    output$infectedvsbystanders_lnc <- renderPlot({
        # generate bins based on input$bins from ui.R
        # draw the histogram with the specified number of bins
        FDRTHRESHOLD <- input$fdr
        FCTHRESHOLD <- log2(input$fc)
        df_filtered <- fcHurdle_all_infectedVSbystanders[fdr<FDRTHRESHOLD & abs(coef)>FCTHRESHOLD,]
        df_plot_ready <- df_filtered %>% group_by(type) %>% dplyr::mutate(count= n_distinct(primerid))
        df_plot_ready <- df_plot_ready[df_plot_ready$type == "lncRNA",]
        ggplot(df_plot_ready, aes(x = type, y =count, fill = type))+
            geom_bar(stat = "identity", color = "black", position =position_dodge())+
            theme_classic()+
            scale_fill_brewer(palette = "Set2")+
            ylim(0,as.integer(max(df_plot_ready$count)))+
            labs(x = "", y = "# DE Genes ", title = "Infected vs Bystanders in Myeloid 24H")
    })
    
    output$celltypes <- renderPlot({
        # generate bins based on input$bins from ui.R
        # draw the histogram with the specified number of bins
        FDRTHRESHOLD <- input$fdr
        FCTHRESHOLD <- log2(input$fc)
        df_filtered <- hurdle_summary_celltypes_lnc[fdr<FDRTHRESHOLD & abs(coef)>FCTHRESHOLD,]
        df_filtered$celltype
        df_plot_ready <- df_filtered %>% group_by(celltype, comparison) %>% mutate(count= n_distinct(primerid))
        ggplot(df_plot_ready, aes(x = celltype, y =count, fill = comparison))+
            geom_bar(stat = "identity", color = "black", position =position_dodge())+
            theme_classic()+
            scale_fill_brewer(palette = "Set2")+
            labs(x = "", y = "# DE Genes ", title = "DE lncRNAs across celltype")
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
