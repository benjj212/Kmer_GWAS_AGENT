library(shiny)
library(ggplot2)
library(gridExtra)
library(plotly)
# Sample function to create a Manhattan plot
create_manhattan_plot <- function(index_file, plot_file, chr, start_pos = NULL, end_pos = NULL, isolate, blast_data=NULL) {
  #renaming the columns names 
  colnames(index_file) <- c("Chromosome", "length")
  colnames(plot_file) <- c("seq","Chromosome", "Position", "pvalue", "color")
  index_file <- index_file[which(index_file[,2] > 500000),]
  xlim_x <- index_file[grep(chr, index_file[,1]),2]
  selected_chromosome <- plot_file[grep(chr, plot_file[,2]),]
  
  if (!is.null(start_pos) && !is.null(end_pos)) {
    selected_chromosome <- selected_chromosome[selected_chromosome$Position >= start_pos & selected_chromosome$Position <= end_pos, ]
  }
  #### setting color labels 
  color_mapping_table <- as.data.frame(cbind(c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#8dd3c7","#bebada","#b3de69","#fccde5","#d9d9d9", "#4d4d4d"),
                               c("A","B", "C","D","E","F", "G", "H", "I", "J", "K", "L", "M", "N", "Z")))
  colnames(color_mapping_table) <- c("color_code", "color_name")
  selected_chromosome <- merge(selected_chromosome, color_mapping_table, by.x="color", by.y="color_code", all.x=T)
  print(head(selected_chromosome))
  print(head(color_mapping_table))
  gg <- ggplot() + 
    geom_point(data= selected_chromosome, aes(x = Position, y = -log10(pvalue), color=color)) +
    scale_color_manual(values = unique(selected_chromosome$color), labels=unique(selected_chromosome$color_name)) +
    labs(title = paste("Manhattan Plot - Chromosome", chr,"-",isolate),
         x = "Position", y = "-log10(P-value)") +
    theme(axis.text.x=element_text(size=15)) +
    xlim(0, xlim_x) +
    ylim(3,25)
  
  if (!is.null(blast_data) && nrow(blast_data) > 0) {
    blast_data$ymin <- as.numeric(blast_data$ymin)  # Convert ymin to numeric if not already
    blast_data$ymax <- as.numeric(blast_data$ymax)  # Convert ymax to numeric if not already
    blast_data$start <- as.numeric(blast_data$start)
    blast_data$end <- as.numeric(blast_data$end)
    if(sum(blast_data$start)>0){
      gg <- gg + geom_rect(data=blast_data, aes(xmin=start, xmax=end, ymin=ymin, ymax=ymax), linewidth=1, color="black") + 
        geom_label(data=blast_data, aes(x=start, y=ymin-1, label=Pmgene))
    }
  }
  if (!is.null(start_pos) && !is.null(end_pos)) {
    gg <- gg + xlim(start_pos, end_pos)
  }
  return(gg)
}

# Sample data for genomes and isolates
genomes <- c("Chinese_Spring_v2.1", "Taes_Fielder_v1_genome", "Taes_Jagger_v1_1_genome", "Taes_Julius_v1_genome", "Taes_Lancer_v1_genome", 
             "Taes_Landmark_v1_v1_genome", "Taes_Mace_v1_genome", "Taes_Renan_v1_genome", "Taes_Stanley_v1_2_genome", "Taes_SYMattis_v1_1_genome", 
             "Genome_Aegilops_taushii_strangulata", "Genome_Triticum_urartu","Triticum_spelta.PGSBv2.0.dna.toplevel", "Triticum_turgidum.Svevo.v1.dna.toplevel")
isolates <- c("ARG_4.2","CHE_96224", "CHE_97251", "GRB_JIW2", "IRN_GOR5", "JPN.Chikara", "KAZ_1b", "POL_3", "THUN12", "TUR_1C")
Chromosomes <- c("1A","1B","1D","2A","2B","2D","3A","3B","3D","4A","4B","4D","5A","5B","5D","6A","6B","6D","7A","7B","7D")

# Define UI
ui <- fluidPage(
  titlePanel("Genome and Isolate"),
  # Main layout with a fluid row
  fluidRow(
    # Sidebar layout with input and output definitions
    column(width = 3,  # Specify width as 3 out of 12 columns, i.e., one-quarter of the width
           sidebarPanel(width=12,
                        # File input for uploading user data
                        fileInput("userfile", "Choose CSV File",
                                  multiple = FALSE,
                                  accept = c("text/csv", 
                                            "text/comma-separated-values,text/plain", 
                                             ".csv")),
                        tags$hr(),
                        # Select input for genomes
                        selectInput("genome",
                                    label = "Select a genome:",
                                    choices = genomes),
                        
                        # Select input for isolates
                        selectInput("isolates",
                                    label = "Select isolates:",
                                    choices = isolates,
                                    multiple = TRUE),
                        
                        # Select input for chromosome
                        selectInput("chromosome",
                                    label = "Select chromosome:",
                                    choices = c("1A","1B","1D","2A","2B","2D","3A","3B","3D","4A","4B","4D","5A","5B","5D","6A","6B","6D","7A","7B","7D")),
                        
                        # Checkbox input for zoom feature

                        textInput("start_pos",
                                  label = "Start Position:",
                                  value = ""),
                        textInput("end_pos",
                                  label = "End Position:",
                                  value = ""),
                        checkboxInput("zoom",
                                      label = "Zoom In",
                                      value = FALSE),
                        
           )
    ),
    
    # Main panel for displaying results
    column(width = 9,  # Specify width as 9 out of 12 columns, i.e., three-quarters of the width
           mainPanel(
             textOutput("selected_genome"),
             textOutput("selected_isolates"),
             plotOutput("manhattan_plots")
           )
    )
  ), 
  fluidRow(
    column(width = 5,
           mainPanel(
             img(src='Upset_plot_app.png', height="100%", width="100%", align="left")
           )
    )
  )
)

server <- function(input, output) {
  
  # Output the selected genome
  output$selected_genome <- renderText({
    paste("Selected genome:", input$genome)
  })
  
  # Output the selected isolates
  output$selected_isolates <- renderText({
    paste("Selected isolates:", paste(input$isolates, collapse = ", "))
  })
  
  # Reactive expression to read uploaded data
  uploaded_data <- reactive({
    req(input$userfile)
    df <- read.table(input$userfile$datapath, header = TRUE, sep = " ")[, c(5, 1, 2, 4)]
    
    # Ensure the dataframe has the required columns for a Manhattan plot
    validate(
      need(all(c("seq","chr", "start", "pval") %in% names(df)), "Data must contain columns: 'seq','chr', 'start', 'pval'")
    )
    
    return(df)
  })
  
  # Render Manhattan plots
  output$manhattan_plots <- renderPlot({
    if (!is.null(input$isolates)){
      plots <- lapply(input$isolates, function(isolate){
        print(input$isolates)
        index_file_path <- paste0(input$genome, ".fasta.fai")
        if (!file.exists(index_file_path)) {
          print(paste("Index file not found:", index_file_path))
          return(NULL)
        }
        index_file <- read.table(index_file_path)[, c(1, 2)]
        blast_matrix <- read.table("matrix.all.txt", header=T, as.is=T)
        blast_matrix_selected <- blast_matrix[grep(input$genome, blast_matrix[,1]),]
        blast_matrix_selected <- blast_matrix_selected[which(as.numeric(blast_matrix_selected$start)>0),]
        blast_matrix_selected <- blast_matrix_selected[grep(input$chromosome, blast_matrix_selected$Chr),]
        if(length(blast_matrix_selected[,1])==0){
          blast_matrix_selected <- matrix(c(input$genome, "nogene","nochrome", 0,0), nrow = 1, ncol = 5)
        }
        print("test_1")
        blast_matrix_selected <- cbind(blast_matrix_selected,4,5)
        print(blast_matrix_selected)
        colnames(blast_matrix_selected) <- c("genome", "Pmgene", "Chr", "start", "end", "ymin", "ymax")
        print("test3")
        print(head(blast_matrix_selected))
        blast_matrix_selected <- as.data.frame(blast_matrix_selected)
        print(blast_matrix_selected)
        GWAS_data_path <- paste0("Ci_2.", isolate, "_", input$genome, ".txt")
        if (!file.exists(GWAS_data_path)) {
          print(paste("GWAS data file not found:", GWAS_data_path))
          return(NULL)
        }
        
        #### reading the colors set up 
        first <- read.table("first_bar.txt", as.is=T)
        second <- read.table("second_bar.txt", as.is=T)
        third <- read.table("third_bar.txt", as.is=T)
        fourth <- read.table("fourth_bar.txt", as.is=T)
        five <- read.table("five_bar.txt", as.is=T)
        six <- read.table("six_bar.txt", as.is=T)
        seven <- read.table("seven_bar.txt", as.is=T)
        eight <- read.table("eight_bar.txt", as.is=T)
        nine <- read.table("nine_bar.txt", as.is=T)
        ten <- read.table("ten_bar.txt", as.is=T)
        eleven <- read.table("eleven_bar.txt", as.is=T)
        twelve <- read.table("twelve_bar.txt", as.is=T)
        thirteen <- read.table("thirteen_bar.txt", as.is=T)
        fourteen <- read.table("fourteen_bar.txt", as.is=T)
        
        # Use the uploaded data if available
        if (!is.null(input$userfile)) {
          GWAS_data <- uploaded_data()
        } else {
          GWAS_data <- read.table(GWAS_data_path)[, c(5, 1, 2, 4)]
        }
        
        GWAS_data$color <- "#4d4d4d"
        GWAS_data[which(GWAS_data[,1] %in% first[,1]),5] <- "#a6cee3"
        GWAS_data[which(GWAS_data[,1] %in% second[,1]),5] <- "#1f78b4"
        GWAS_data[which(GWAS_data[,1] %in% third[,1]),5] <- "#b2df8a"
        GWAS_data[which(GWAS_data[,1] %in% fourth[,1]),5] <- "#33a02c"
        GWAS_data[which(GWAS_data[,1] %in% five[,1]),5] <- "#fb9a99"
        GWAS_data[which(GWAS_data[,1] %in% six[,1]),5] <- "#e31a1c"
        GWAS_data[which(GWAS_data[,1] %in% seven[,1]),5] <- "#fdbf6f"
        GWAS_data[which(GWAS_data[,1] %in% eight[,1]),5] <- "#ff7f00"
        GWAS_data[which(GWAS_data[,1] %in% nine[,1]),5] <- "#cab2d6"
        GWAS_data[which(GWAS_data[,1] %in% ten[,1]),5] <- "#8dd3c7"
        GWAS_data[which(GWAS_data[,1] %in% eleven[,1]),5] <- "#bebada"
        GWAS_data[which(GWAS_data[,1] %in% twelve[,1]),5] <- "#b3de69"
        GWAS_data[which(GWAS_data[,1] %in% thirteen[,1]),5] <- "#fccde5"
        GWAS_data[which(GWAS_data[,1] %in% fourteen[,1]),5] <- "#d9d9d9"
        
        
        # Check if zoom feature is activated
        if (input$zoom) {
          # Check if start and end positions are provided
          if (!is.null(input$start_pos) && !is.null(input$end_pos)) {
            create_manhattan_plot(index_file, GWAS_data, input$chromosome, as.numeric(input$start_pos), as.numeric(input$end_pos), isolate=isolate, blast_data=blast_matrix_selected)
          } else {
            # If start and end positions are not provided, show entire genome
            create_manhattan_plot(index_file, GWAS_data, input$chromosome, isolate=isolate, blast_data=blast_matrix_selected)
          }
        } else {
          # If zoom feature is not activated, show entire genome
          create_manhattan_plot(index_file, GWAS_data, input$chromosome, isolate=isolate, blast_data=blast_matrix_selected)
        }
      })
      if (length(plots) == 0) {
        plot(NULL, xlim = c(0, 1), ylim = c(0, 1), main = "No plots to display")
      } else {
        do.call(grid.arrange, c(plots, ncol = 1))
      }
    } else {
      plot(NULL, xlim = c(0, 1), ylim = c(0, 1), main = "No isolates selected")
    }
  }, height=900, width=1400)
}
# Run the application
shinyApp(ui = ui, server = server)