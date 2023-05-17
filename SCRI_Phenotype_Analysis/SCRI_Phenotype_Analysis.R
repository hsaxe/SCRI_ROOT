#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(shinythemes)
library(shinydashboard)
library(shinyjs)
library(tidytext)
library(data.table)
library(DT)
library(stringr)
library(shinyWidgets)
library(ggpubr)

dat = read.csv('SCRI_Phenotype_Data_InComplete_Cases.csv')

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("SCRI Phenotype Analysis"),
  
  # Sidebar with a select inputs for multiple variables
  sidebarLayout(
    sidebarPanel(
      
      selectInput('plotType',
                  'Select plot type',
                  choices = NULL),
      
      selectInput('x',
                  'Select x-axis value/factor',
                  choices = NULL),
      
      selectInput('y',
                  'Select y-axis value/factor',
                  choices = NULL),
      
      checkboxInput('transform',
                    'Log10 transform nematode data'),
      
      # selectInput('fill',
      #             'Select what to fill color by',
      #             choices = NULL),
      # 
      # selectInput('facet1',
      #             'Select facet one',
      #             choices = NULL),
      # 
      # selectInput('facet2',
      #             'Select facet two',
      #             choices = NULL),
      # 
      # pickerInput('subset',
      #             'Subset x-axis data',
      #             options = list(`actions-box` = TRUE),
      #             choices = NULL,
      #             multiple = T),
      # 
      # selectInput('filterCol',
      #             'Filter dataset: Select column',
      #             choices = NULL,
      #             multiple = F),
      # 
      # pickerInput('filterVal',
      #             'Filter dataset: Select value(s)',
      #             options = list(`actions-box` = TRUE),
      #             choices = NULL,
      #             multiple = T),
      # 
      # selectInput('filterCol2',
      #             'Filter dataset again: Select column',
      #             choices = NULL,
      #             multiple = F),
      # 
      # pickerInput('filterVal2',
      #             'Filter dataset again: Select value(s)',
      #             options = list(`actions-box` = TRUE),
      #             choices = NULL,
      #             multiple = T),
      
      textInput('Plot_name', 'Name this plot'),
      
      numericInput('Save_width', 'Change size of saved plot',
                   value = 8),
      
      downloadButton('save_plot', 
                     'Save this plot'),
      
      width = 2
    ),
    
    # Show a plot of the selected variables and data table used to generate it
    mainPanel(
      plotOutput("plot",
                 height = '600px'),
      
      br(),
      
      DT::DTOutput('plotData'),
      
      width = 10
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  # Dataset must be reactive to accept inputs from ui
  dataset = reactive({
    
    dat
    
  })
  
  # datasetMod = reactiveValues()
  # 
  # observe({
  # 
  #   datasetMod$x = dataset()
  # 
  # })
  
  # Renders data table to ui
  output$plotData = DT::renderDT({
    
    datatable(dataset2(), 
              selection = list(target = 'column'),
              filter = 'top')
  })
  
  
  # These observers are needed for dynamic plotting of variables
  
  observe({
    updateSelectInput(session, 'plotType',
                      'Select plot type',
                      choices = c('Scatter', 'Density'),
                      selected = 'Scatter')
  })
  
  observe({
    updateSelectInput(session, 'x',
                      'Select x-axis value/factor',
                      choices = names(dataset()),
                      selected = 'CG_Avg')
  })
  
  observe({
    updateSelectInput(session, 'y',
                      'Select y-axis value/factor',
                      choices = names(dataset()),
                      selected = 'X2Y_length')
  })
  
  dataset2 = reactive({
    
    if(input$transform == T){
      
      dataset = dataset() %>% 
        mutate(across(matches('RLN'), ~ log10(.x+1)))
      
    } else {
      
      dataset()
      
    } 
    
  })
  

  
  # Making reactive plot
  plotInput = reactive({
    
    if(input$plotType == 'Scatter') {
      
      ggplot(dataset2(),
             aes(.data[[input$x]], .data[[input$y]],
                 fill = factor(Female_Parent),
                 color = factor(Female_Parent)))+
        geom_point(size = 3)+
        geom_smooth(method = 'lm')+
        # stat_regline_equation()+
        # stat_cor(method = 'spearman',
        #          size = 6)+
        stat_smooth(method = 'lm')+
        # stat_regline_equation(label.y = dataset() %>% 
        #                         pull(.data[[input$y]]) %>% 
        #                         na.omit() %>% 
        #                         max() * 0.95,
        #                       aes(label =  paste(after_stat(eq.label), ..adj.rr.label.., sep = "~~~~")),
        #                       size = 6)+
        stat_cor(method = 'spearman', 
                 size = 6)+
        theme_gray(base_size = 15)
      
    } else {
      
      ggplot(dataset2(),
             aes(.data[[input$x]],
                 color = factor(Female_Parent)))+
        geom_density(size = 2)+
        theme_gray(base_size = 15)
      
    }
    
    
    
    
    
  })
  
  # Rendering reactive plot to ui
  output$plot = renderPlot({
    
    print(plotInput())
    
  })
  
  ## This creates a button to downnload the current plot
  output$save_plot <- downloadHandler(
    filename = function() { paste(input$Plot_name, '.png', sep='_') },
    content = function(file) {
      save_plot(file, plotInput(), base_height = input$Save_width)
    }
  )
}

# Run the application 
shinyApp(ui = ui, server = server)
