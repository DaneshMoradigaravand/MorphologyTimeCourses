require(shiny)
library(tidyverse)

sigmoid = function(params, x) {
  params[1] / (1 + exp(-params[2] * (x - params[3])))
}

file_tags<-read_csv("gene_locus_file.csv", show_col_types = FALSE)
time_points_files<-read_csv("PAmorphology.csv", show_col_types = FALSE)
time_points_files_colour<-read_csv("PAcolour.csv", show_col_types = FALSE)

annotation<-read_csv("annotation.csv", show_col_types = FALSE)
properties_morphology<-read.csv("Properties_morphology.csv")
properties_colony_size<-read.csv("Properties_colony_size.csv")
properties_colony_colour<-read.csv("Properties_colony_colour.csv")

plot_gene_from_file_morophology<-function(input_gene,input,properties){
  unique_loci<-unique(input$PA14.ID)
  timepoints<-c("13h", "19h", "39h", "44h", "69h", "89h")
  x_timepoints<-match(input$timepoint[input$PA14.ID==input_gene],timepoints )
  x_timepoints<-x_timepoints+runif(n=length(x_timepoints), min=1e-12, max=.01)
  y_growth<-input$morphology.score.fixed.circles[input$PA14.ID==input_gene]
  
  params<-properties[properties$locus.id==input_gene,c(2,3,4)]
  as.numeric(params)
  
  x_timepoint_range<-seq(min(x_timepoints),max(x_timepoints),0.1)
  y_growth_predicted <- sigmoid(as.numeric(params),x_timepoint_range)
  
  plot(x_timepoints,y_growth)
  df<-data.frame(list(x=x_timepoints, y=y_growth))
  df_predict<-data.frame(list(x=x_timepoint_range, y=y_growth_predicted ))
  
  if(   !is.na(as.numeric(params)[1]) ){
  output<-ggplot(df,aes(x,y))+
    geom_point(col="red")+
    geom_line(data=df_predict,size=2)+ 
    ggtitle(paste0("Gene ", input_gene))+
    theme_bw()+
    ylim(range(0,500))+
    ylab("Morphology")+
    xlab("Time Pionts")+ 
    theme(axis.text.x = element_text( size=11, angle=45,hjust = 1),
          axis.text.y = element_text( size=11, hjust = 1),
          axis.title.x = element_text(color="black", size=14, face="bold"),
          axis.title.y = element_text(color="black", size=14, face="bold"),
          strip.text.x = element_text(
            size = 13, color = "black", face = "bold"
          ))+
    scale_x_continuous(breaks=1:6,labels=c("13h", "19h", "39h","44h", "69h", "89h"))
  return(output)
  }
  else{
    output<-ggplot(df,aes(x,y))+
      geom_point(col="red")+
      ggtitle(paste0("Gene ", input_gene))+
      theme_bw()+
      ylim(range(0,500))+
      ylab("Morphology")+
      xlab("Time Pionts")+ 
      theme(axis.text.x = element_text( size=11, angle=45,hjust = 1),
            axis.text.y = element_text( size=11, hjust = 1),
            axis.title.x = element_text(color="black", size=14, face="bold"),
            axis.title.y = element_text(color="black", size=14, face="bold"),
            strip.text.x = element_text(
              size = 13, color = "black", face = "bold"
            ))+
      scale_x_continuous(breaks=1:6,labels=c("13h", "19h", "39h","44h", "69h", "89h"))
    return(output)
    
  }
  
  
}
plot_gene_from_file_colony_size<-function(input_gene,input,properties){
  unique_loci<-unique(input$PA14.ID)
  timepoints<-c("13h", "19h", "39h", "44h", "69h", "89h")
  x_timepoints<-match(input$timepoint[input$PA14.ID==input_gene],timepoints )
  x_timepoints<-x_timepoints+runif(n=length(x_timepoints), min=1e-12, max=.01)
  y_growth<-input$colony.size[input$PA14.ID==input_gene]
  
  params<-properties[properties$locus.id==input_gene,c(2,3,4)]
  as.numeric(params)
  
  x_timepoint_range<-seq(min(x_timepoints),max(x_timepoints),0.1)
  y_growth_predicted <- sigmoid(as.numeric(params),x_timepoint_range)
  
  plot(x_timepoints,y_growth)
  df<-data.frame(list(x=x_timepoints, y=y_growth))
  df_predict<-data.frame(list(x=x_timepoint_range, y=y_growth_predicted ))
  
  if(   !is.na(as.numeric(params)[1]) ){
    output<-ggplot(df,aes(x,y))+
      geom_point(col="red")+
      geom_line(data=df_predict,size=2)+ 
      ggtitle(paste0("Gene ", input_gene))+
      theme_bw()+
      ylim(range(0,50000))+
      ylab("Colony Size")+
      xlab("Time Pionts")+ 
      theme(axis.text.x = element_text( size=11, angle=45,hjust = 1),
            axis.text.y = element_text( size=11, hjust = 1),
            axis.title.x = element_text(color="black", size=14, face="bold"),
            axis.title.y = element_text(color="black", size=14, face="bold"),
            strip.text.x = element_text(
              size = 13, color = "black", face = "bold"
            ))+
      scale_x_continuous(breaks=1:6,labels=c("13h", "19h", "39h","44h", "69h", "89h"))
    return(output)
  }
  else{
    output<-ggplot(df,aes(x,y))+
      geom_point(col="red")+
      ggtitle(paste0("Gene ", input_gene))+
      theme_bw()+
      ylim(range(0,50000))+
      ylab("Colony Size")+
      xlab("Time Pionts")+ 
      theme(axis.text.x = element_text( size=11, angle=45,hjust = 1),
            axis.text.y = element_text( size=11, hjust = 1),
            axis.title.x = element_text(color="black", size=14, face="bold"),
            axis.title.y = element_text(color="black", size=14, face="bold"),
            strip.text.x = element_text(
              size = 13, color = "black", face = "bold"
            ))+
      scale_x_continuous(breaks=1:6,labels=c("13h", "19h", "39h","44h", "69h", "89h"))
    return(output)
    
  }
  
  
}
plot_gene_from_file_colour<-function(input_gene,input,properties){
  unique_loci<-unique(input$PA14.ID)
  timepoints<-c("13h", "19h", "39h", "44h", "69h", "89h")
  x_timepoints<-match(input$time_point[input$gene==input_gene],timepoints )
  x_timepoints<-x_timepoints+runif(n=length(x_timepoints), min=1e-12, max=.01)
  y_growth<-input$value[input$gene==input_gene]
  
  params<-properties[properties$locus.id==input_gene,c(2,3,4)]
  as.numeric(params)
  
  x_timepoint_range<-seq(min(x_timepoints),max(x_timepoints),0.1)
  y_growth_predicted <- sigmoid(as.numeric(params),x_timepoint_range)
  
  plot(x_timepoints,y_growth)
  df<-data.frame(list(x=x_timepoints, y=y_growth))
  df_predict<-data.frame(list(x=x_timepoint_range, y=y_growth_predicted ))
  
  if(   !is.na(as.numeric(params)[1]) ){
    output<-ggplot(df,aes(x,y))+
      geom_point(col="red")+
      geom_line(data=df_predict,size=2)+ 
      ggtitle(paste0("Gene ", input_gene))+
      theme_bw()+
      ylim(range(0,20000000))+
      ylab("Colony Colour")+
      xlab("Time Pionts")+ 
      theme(axis.text.x = element_text( size=11, angle=45,hjust = 1),
            axis.text.y = element_text( size=11, hjust = 1),
            axis.title.x = element_text(color="black", size=14, face="bold"),
            axis.title.y = element_text(color="black", size=14, face="bold"),
            strip.text.x = element_text(
              size = 13, color = "black", face = "bold"
            ))+
      scale_x_continuous(breaks=1:6,labels=c("13h", "19h", "39h","44h", "69h", "89h"))
    return(output)
  }
  else{
    output<-ggplot(df,aes(x,y))+
      geom_point(col="red")+
      ggtitle(paste0("Gene ", input_gene))+
      theme_bw()+
      ylim(range(0,20000000))+
      ylab("Colony Colour")+
      xlab("Time Pionts")+ 
      theme(axis.text.x = element_text( size=11, angle=45,hjust = 1),
            axis.text.y = element_text( size=11, hjust = 1),
            axis.title.x = element_text(color="black", size=14, face="bold"),
            axis.title.y = element_text(color="black", size=14, face="bold"),
            strip.text.x = element_text(
              size = 13, color = "black", face = "bold"
            ))+
      scale_x_continuous(breaks=1:6,labels=c("13h", "19h", "39h","44h", "69h", "89h"))
    return(output)
    
  }
  
  
}

ranking_morphology<-read.csv("property_morphology_tot.csv")  %>% mutate(feature="Morphology")
ranking_colony_size<-read.csv("property_colony_size_tot.csv")  %>% mutate(feature="Colony Size")
ranking_colony_colour<-read.csv("property_colony_colour_tot.csv")  %>% mutate(feature="Colony Colour")


ui <- fluidPage(
    headerPanel("Morphology Screen"),
    sidebarPanel(
        helpText("Select a gene locus, the output will be the images of clones, gene information and measured attributes"),
        selectInput("gene", 
                    label = "Choose a Gene",
                   # choices = NULL,
                    choices = gsub("\\.jpg","", file_tags$Tags),
                   selected = "PA14_73390")
    ),
    mainPanel(
              tableOutput('to'),
              tableOutput('col'),
              fluidRow(
                column(4,plotOutput('plot_morph')), 
                column(4,plotOutput('plot_size')),
                column(4,plotOutput('plot_colour'))
              ),
              fluidRow(
              
              column(1, offset = 3,plotOutput("plot1"))
              )
              #column(4, offset = 1,)
              
              
              
    )
)
server <- function(input, output) {

    output$to <- renderTable(annotation[match(input$gene,annotation$Gene_Locus),])
    #output$morph <- renderTable(ranking_morphology[match(input$gene,ranking_morphology$Gene),] )
    #output$size <- renderTable(ranking_colony_size[match(input$gene,ranking_colony_size$Gene),] )
    output$col <- renderTable(
      rbind( 
        ranking_morphology[match(input$gene,ranking_morphology$Gene),] ,
        ranking_colony_size[match(input$gene,ranking_colony_size$Gene),],
        ranking_colony_colour[match(input$gene,ranking_colony_colour$Gene),])
      )
    
    
    output$plot1 <- renderImage({ 
      filename <- normalizePath(file.path('./www',paste0(input$gene,".jpg")))
      
      list(src = filename,
           contentType = 'image/jpg',
           alt = "This is alternate text")
    }, deleteFile = FALSE)
    
    output$plot_morph <- renderPlot({
      plot_gene_from_file_morophology(input$gene, time_points_files, properties_morphology)
    })
    
    output$plot_size <- renderPlot({
      plot_gene_from_file_colony_size(input$gene, time_points_files, properties_colony_size)
    })
    
    output$plot_colour <- renderPlot({
      plot_gene_from_file_colour(input$gene, time_points_files_colour, properties_colony_colour)
    })
    
}

shinyApp(ui, server, options = list(height = 1080))

