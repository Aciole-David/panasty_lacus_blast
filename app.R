library(shiny)
require(XML)
library(plyr)
library(dplyr)
library(shinythemes)
library(DT)
library(sfsmisc)
library(Biostrings)
library(shinyjs)

# Define ports to expose
options(shiny.host = "0.0.0.0")
options(shiny.port = 8080)

pcinfo<-as.list(Sys.cpuinfo())
threads<-(as.numeric(pcinfo$siblings)-2)

#ref<-readDNAStringSet('astyanax-pan.fasta')
#rnm<-data.frame(headers=names(ref))


# system("wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-linux.tar.gz")
# system("tar xvpf ncbi-blast-2.16.0+-x64-linux.tar.gz")
# system("cp ncbi-blast-2.16.0+/bin/blastn .")
# system("cp ncbi-blast-2.16.0+/bin/tblastn .")
system("chmod +x /home/shiny-app/www/blastn")
system("chmod +x /home/shiny-app/www/tblastn")
system("chmod +x /home/shiny-app/www/seqkit")

custom_db <- c("Pan Astyanax lacustris")


ui <- fluidPage(
  #rclipboardSetup(),
  theme = shinytheme("cerulean"),
                tagList(
                  tags$head(
                    tags$link(rel="stylesheet", type="text/css",href="style.css"),
                    tags$script(type="text/javascript", src = "busy.js")
                  )
                ),
                
                #This block gives us all the inputs:
                mainPanel(width = 12,
                  #uiOutput("clip"),
                  h4('Astyanax lacustrix Pantranscriptome BLAST'),
                  textAreaInput('query', 'Input query sequence below:', value = "", placeholder = "", width = "100%", height="100px"),
                  div(style="display:inline-block",
                  selectInput("db", "Database:", choices=c(custom_db,"nr"), width="190px")),
                  div(style="display:inline-block",
                      selectInput("program", "Program:", choices=c("blastn","tblastn"), width="100px")),
                  #task
                  div(style="display:inline-block",
                      selectInput("task", "Task:",selected = 'blastn',
                                  choices=c('blastn','blastn-short','dc-megablast','megablast','rmblastn',
                                            'tblastn', 'tblastn-fast'),
                                  width="130px")),
                  div(style="display:inline-block",
                      selectInput("eval", "e-value:", choices=c(1,0.001,1e-4,1e-5,1e-10), width="120px")),
                  actionButton("blast", "BLAST!")
                ),
                
                #this snippet generates a progress indicator for long BLASTs
                #div(class = "busy",  
                #    p("Calculation in progress.."), 
                #    #tags$img(src="https://i.stack.imgur.com/8puiO.gif", height = 100, width = 100,align = "center")
                #    tags$img(src="8puiO.gif", height = 100, width = 100,align = "center"),
                #  ),
                
                #Basic results output
                mainPanel(width = 12,
                  h4("Results"),
                  DT::dataTableOutput("blastResults"),
                  p("",
                    #tableOutput("clicked")
                    ),
                  p('Alignment'),
                  verbatimTextOutput("alignment"),
                  p('Hit sequence'),
                  div(style="height:50px",
                  verbatimTextOutput("alignment2")
                  )
                  )
)

server <- function(input, output, session){
  shinyjs::hide(id = "clicked")

  custom_db <- c("Pan Astyanax lacustris")
  custom_db_path <- c("/home/shiny-app/www/astypan")
  
  #task <- input$task


  blastresults <- eventReactive(input$blast, {
    
    #gather input and set up temp file
    query <- input$query
    tmp <- tempfile(fileext = ".fa")
    #task <- input$task

    
    #if else chooses the right database
    if (input$db == custom_db){
      db <- custom_db_path
      remote <- c("")
    } else {
      db <- c("nr")
      #add remote option for nr since we don't have a local copy
      remote <- c("-remote")
    }
    
    
    #this makes sure the fasta is formatted properly
    if (startsWith(query, ">")){
      writeLines(query, tmp)
    } else {
      writeLines(paste0(">Query\n",query), tmp)
    }
    
    #calls the blast
    # system(paste0("./",input$program," -query ",tmp," -db ",db," -evalue ",input$eval,
    #               " -task ",input$task, ' -num_threads ', threads,
    #               " -outfmt 5 -max_hsps 1 -max_target_seqs 10 ",remote,'> blast.xml'), intern = T)
    # 
    data <- system(paste0("/home/shiny-app/www/",input$program," -query ",tmp," -db ",db," -evalue ",input$eval,
                          " -task ",input$task, ' -num_threads ', threads, 
                          " -outfmt 5 -max_hsps 1 -max_target_seqs 10 ",remote), intern = T)
    
    xmlParse(data)
    }, ignoreNULL= T)
  

  #Now to parse the results...
  parsedresults <- reactive({
    if (is.null(blastresults())){}
    else {
      xmltop = xmlRoot(blastresults())
      
      #the first chunk is for multi-fastas
      results <- xpathApply(blastresults(), '//Iteration',function(row){
        query_ID <- getNodeSet(row, 'Iteration_query-def') %>% sapply(., xmlValue)
        hit_IDs <- getNodeSet(row, 'Iteration_hits//Hit//Hit_id') %>% sapply(., xmlValue)
        hit_Def <- getNodeSet(row, 'Iteration_hits//Hit//Hit_def') %>% sapply(., xmlValue)
        hit_length <- getNodeSet(row, 'Iteration_hits//Hit//Hit_len') %>% sapply(., xmlValue)
        bitscore <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_bit-score') %>% sapply(., xmlValue)
        eval <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_evalue') %>% sapply(., xmlValue)
        cbind(query_ID,hit_IDs,hit_Def,hit_length,bitscore,eval)
      })
      
      #this ensures that NAs get added for no hits
      results <-  rbind.fill(lapply(results,function(y){as.data.frame((y),stringsAsFactors=FALSE)}))
      results2 <<- results
    }
  })
  
  #makes the datatable
  output$blastResults <- renderDataTable(caption = "Select a hit to view alignments",{
    if (is.null(blastresults())){}
    #better aesthetics table and help when no hits
    else if (length(parsedresults()) < 2 ){
      data.frame('No Hits found'='Try changing parameters',
                 row.names = NULL, check.names = F)
    }
    else {
      #print(str(parsedresults()))
      #print(length(parsedresults()))
      #print('parsed results ELSE')
      parsedresults()
    }
    
  } , selection="single")
  
  #this chunk gets the alignemnt information from a clicked row
  output$clicked <- renderTable(caption = "", caption.placement = "top",{
    
    if(is.null(input$blastResults_rows_selected)){}
    else if (length(parsedresults()) < 2 ) {
    }
    else{
      xmltop = xmlRoot(blastresults())
      clicked = input$blastResults_rows_selected
      tableout<- data.frame(parsedresults()[clicked,])
      
      tableout <- t(tableout)
      names(tableout) <- c("")
      rownames(tableout) <- c("Query ID","Hit ID","Hit Def", "Length", "Bit Score", "e-value")
      colnames(tableout) <- c("Alignment\ selected")
      data.frame(tableout, check.names = F)
    }
  },rownames =T,colnames =T)
  
  #this chunk makes the alignments for clicked rows
  output$alignment <- renderText({
    if(is.null(input$blastResults_rows_selected)){}
    else if (length(parsedresults()) < 2 ) {
      #'No hits found'
    }
    else {
      xmltop = xmlRoot(blastresults())
      
      clicked = input$blastResults_rows_selected

      #loop over the xml to get the alignments
      align <<- xpathApply(blastresults(), '//Iteration',function(row){
        top <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_qseq') %>% sapply(., xmlValue)
        mid <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_midline') %>% sapply(., xmlValue)
        bottom <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_hseq') %>% sapply(., xmlValue)
        rbind(top,mid,bottom)
      })

      #split the alignments every 40 carachters to get a "wrapped look"
      alignx <<- do.call("cbind", align)
      splits <- strsplit(gsub("(.{80})", "\\1,", alignx[1:3,clicked]),",")

      #paste them together with returns '\n' on the breaks
      split_out <- lapply(1:length(splits[[1]]),function(i){
        rbind("\n",
              paste0("Q - ",splits[[1]][i],"\n"),
              paste0("M - ",splits[[2]][i],"\n"),
              paste0("H - ",splits[[3]][i],"\n"))
      })
      unlist(split_out)
    }
  })
  
  # // Get subject ID
    output$clicked2 <- renderTable(caption = "", caption.placement = "top",{
    
    if(is.null(input$blastResults_rows_selected)){}
    else if (length(parsedresults()) < 2 ) {
    }
    else{
      xmltop = xmlRoot(blastresults())
      clicked2 = input$blastResults_rows_selected
      tableout<- data.frame(parsedresults()[clicked2,2])
      
      #tableout <- t(tableout)
      #names(tableout) <- c("")
      #rownames(tableout) <- c("Query ID","Hit ID","Hit Def", "Length", "Bit Score", "e-value")
      #colnames(tableout) <- c("Alignment\ selected")
      data.frame(tableout, check.names = F)
    }
  },rownames =T,colnames =T)
  

    # // REACTIVE GET HIT SEQ
    getseq <- reactive({
       if(is.null(input$blastResults_rows_selected)){}
    else if (length(parsedresults()) < 2 ) {
      #'No hits found'
    }
    else {
      xmltop = xmlRoot(blastresults())

      clicked2 = input$blastResults_rows_selected

      #loop over the xml to get the alignments
      align <<- xpathApply(blastresults(), '//Iteration',function(row){
        top2 <- getNodeSet(row, 'Iteration_hits//Hit//Hit_id') %>% sapply(., xmlValue)
        rbind(top2)
      })

      #split the alignments every 40 carachters to get a "wrapped look"
      alignx <<- do.call("cbind", align)
      splits <- strsplit(gsub("(.{80})", "\\1,", alignx[1:1,clicked2]),",")

      #paste them together with returns '\n' on the breaks
      split_out <- lapply(1:length(splits[[1]]),function(i){
        rbind(splits[[1]][i])
      })
      seqid<-as.character(unlist(split_out))
      print(seqid)
      #contig312.1
      data2<-as.character(system(paste0(
      "/home/shiny-app/www/seqkit grep -w 0 -p \'",seqid,"\' /home/shiny-app/www/astyanax-pan.fasta -j ",threads),
                     intern = T))
      #system(paste0("./seqkit grep -p \'",seqid,"\' astyanax-pan.fasta -j 6"),
      #               intern = T)
      # data2<-system(paste0(
      #   "./seqkit grep -w 0 -p 'contig0.1' hasty.fas -j 6 "),
      #                 intern = T)
      data2h<-data2[1]
      data2s<-data2[2]
      data2ss<-as.character(strsplit(gsub("(.{100})", '\\1\n', data2s), ','))
      paste0(data2h,'\n',data2ss)
      #rbind(hdata2, sdata2)
      
      #longtext='ATCGATGCATGCTAGCATGCTAGCTAGCTAGCTAGCTAGCTAGCTGCTAGCTAG'
      #strsplit(gsub("(.{8})", '\\1\n', longtext), ',')
      #as.character(strsplit(gsub("(.{80})", '\\1\n', data2), ','))
      #write(a, file = 'aaaa.txt')
      
    }

  } ) #, ignoreNULL = T)
    
  # //RENDER SUBJECT SEQ
  output$alignment2 <- renderText({
    if(is.null(input$blastResults_rows_selected)){}
    else if (length(parsedresults()) < 2 ) {
      #'No hits found'
    }
    else {
      getseq()

    }
  })
  
    # // REACTIVE GET ALIGN2 CLIPBOARD
  #'   output$clip <- renderUI({
  #'   if(is.null(input$blastResults_rows_selected)){}
  #'   else if (length(parsedresults()) < 2 ) {
  #'     #'No hits found'
  #'   }
  #'   else {
  #'     clipalign2()
  #' 
  #'   }
  #' })
  
    # Add clipboard buttons
  # clipalign2 <- reactive({
  #   rclipButton(
  #     inputId = "clipbuton",
  #     label = "rclipButton Copy",
  #     clipText = output$clicked,
  #     icon = icon("clipboard"),
  #     tooltip = "Click me... I dare you!",
  #     placement = "top",
  #     options = list(delay = list(show = 2800, hide = 100), trigger = "hover")
  #   )
  # })
  

}

# Run the application 
shinyApp(ui = ui, server = server)

