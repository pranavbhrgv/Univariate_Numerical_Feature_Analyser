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
library(purrr)
library(data.table)
library(car)
library(ggplot2)
library(futile.logger)

# Functions used in the server file
p_value_distr_num <-
  function(final_df,
           column_name,
           ITERATIONS,
           TARGET_POSITIVE,
           TARGET_NEGATIVE,
           SAMPLE_SIZE,
           P_VALUE_LOWER,
           P_VALUE_UPPER) {
    p_values = c()
    flog.info("Step 1")
    for (y in c(1:ITERATIONS)) {
      s = sample(final_df[final_df$target == TARGET_POSITIVE, c(column_name)], SAMPLE_SIZE)
      t = sample(final_df[final_df$target == TARGET_NEGATIVE, c(column_name)], SAMPLE_SIZE)
      temp = t.test(s, t)
      p_values[y] = temp$p.value
    }
    
    sig = sum(p_values <= P_VALUE_LOWER)
    not_sig = sum(p_values >= P_VALUE_UPPER)
    
    return(c(sig, not_sig))
    
  }

#Wilcox Signed Test
p_value_distr_graph_nnorm_num <-
  function(final_df,
           column_name,
           ITERATIONS,
           TARGET_POSITIVE,
           TARGET_NEGATIVE,
           SAMPLE_SIZE,
           P_VALUE_LOWER,
           P_VALUE_UPPER) {
    p_values = c()
    flog.info("Step 1")
    for (y in c(1:ITERATIONS)) {
      s = sample(final_df[final_df$target == TARGET_POSITIVE, c(column_name)], SAMPLE_SIZE)
      t = sample(final_df[final_df$target == TARGET_NEGATIVE, c(column_name)], SAMPLE_SIZE)
      temp = wilcox.test(s, t)
      p_values[y] = temp$p.value
    }
    
    sig = sum(p_values <= P_VALUE_LOWER)
    not_sig = sum(p_values >= P_VALUE_UPPER)
    
    return(c(sig, not_sig))
  }

#Check NA
check_na <- function(vec) {
  s = sum(is.na(vec)) / length(vec) * 100
  return(as.numeric(s[1]))
}

#Single Function called in Server File
fun <-
  function(input,
           MISSINGNESS_PERCENTAGE,
           ITERATIONS,
           TARGET_POSITIVE,
           TARGET_NEGATIVE,
           SAMPLE_SIZE,
           P_VALUE_LOWER,
           P_VALUE_UPPER) {
    feature_list <- colnames(input)
    feature_list <-
      feature_list[!feature_list %in% c("personId", "target")]
    
    feature_completeness = input %>% select(feature_list) %>% map_dfr( ~ check_na(.))
    
    f = data.frame(feature_completeness %>% t())
    f$column_names = row.names(f)
    usable_features = f %>% filter(feature_completeness.....t.. <= MISSINGNESS_PERCENTAGE)
    
    #Parametric Way
    temp = input %>% select(usable_features$column_names) %>% select_if(is.numeric) %>% imap_dfr(function(x, name) {
      flog.info(name)
      p_value_distr_num(
        input,
        name,
        ITERATIONS,
        TARGET_POSITIVE,
        TARGET_NEGATIVE,
        SAMPLE_SIZE,
        P_VALUE_LOWER,
        P_VALUE_UPPER
      )
    })
    
    temp2 = temp %>% t() %>% data.frame()
    
    temp2$column_names = row.names(temp2)
    
    colnames(temp2) <- c("=<0.04 count", "=>0.06 count", "column_names")
    row.names(temp2) <- NULL
    
    table(temp2$`=<0.04 count` > temp2$`=>0.06 count`)
    
    temp2$useful_features = ifelse(temp2$`=<0.04 count` > temp2$`=>0.06 count`, 1, 0)
    temp2$is_numeric = 1
    
    feature_usefulness_type = temp2 %>% full_join(usable_features, by = "column_names")
    
    #Non Parametric Way
    temp = input %>% select(usable_features$column_names) %>% select_if(is.numeric) %>% imap_dfr(function(x, name) {
      flog.info(name)
      p_value_distr_graph_nnorm_num(
        input,
        name,
        ITERATIONS,
        TARGET_POSITIVE,
        TARGET_NEGATIVE,
        SAMPLE_SIZE,
        P_VALUE_LOWER,
        P_VALUE_UPPER
      )
    })
    
    temp2 = temp %>% t() %>% data.frame()
    
    
    temp2$column_names = row.names(temp2)
    
    colnames(temp2) <- c("=<0.04 count", "=>0.06 count", "column_names")
    row.names(temp2) <- NULL
    
    table(temp2$`=<0.04 count` > temp2$`=>0.06 count`)
    
    temp2$useful_features = ifelse(temp2$`=<0.04 count` > temp2$`=>0.06 count`, 1, 0)
    temp2$is_numeric = 1
    
    feature_usefulness_type_np = temp2 %>% full_join(usable_features, by =
                                                       "column_names")
    
    
    #Normality
    temp = input %>% select(feature_usefulness_type$column_names) %>% select_if(is.numeric) %>% imap_dfr(function(x, name) {
      out <- tryCatch({
        x = sample(x, size = min(4999, nrow(input)), replace = FALSE)
        s = shapiro.test(x)
        return(c(name, s$p.value))
      },
      error = function(cond) {
        flog.error(paste("Column caused an error:", name))
        #flog.info("Here's the original warning flog.info:")
        #flog.info(cond)
        # Choose a return value in case of error
        return(c(name, NA))
      },
      warning = function(cond) {
        flog.warn(paste("Column caused a warning:", name))
        # Choose a return value in case of warning
        return(c(name, NULL))
      },
      finally = {
        flog.info(paste("Processed column:", name))
      })
      return(out)
    })
    
    shapiro_test = temp %>% t() %>% data.frame(row.names = NULL, stringsAsFactors = FALSE) %>% rename(column_names =
                                                                                                        X1, p_values_shapiro_test = X2) %>% mutate(p_values_shapiro_test = as.numeric(p_values_shapiro_test))
    
    #Homogeneity
    library(car)
    
    tar = as.factor(input$target)
    temp = input %>% select(usable_features$column_names) %>% select_if(is.numeric) %>% imap_dfr(function(x, name) {
      s = leveneTest(x ~ tar)
      return(c(name, s$`Pr(>F)`[1]))
    })
    
    lev_test = temp %>% t() %>% data.frame(row.names = NULL, stringsAsFactors = FALSE) %>% rename(column_names =
                                                                                                    X1, p_values_levene_test = X2) %>% mutate(p_values_levene_test = as.numeric(p_values_levene_test))
    
    #Merging for all numeric test
    feature_usefulness_type_np = feature_usefulness_type_np %>% rename(useful_features_np =
                                                                         useful_features,
                                                                       column_names_np =
                                                                         column_names)
    feature_usefulness_type_np$feature_completeness.....t.. = NULL
    feature_usefulness_type_np$is_numeric = NULL
    merged_np_p = feature_usefulness_type %>% inner_join(feature_usefulness_type_np,
                                                         by = c("column_names" = "column_names_np"))
    all_merged = merged_np_p %>% full_join(lev_test, by = "column_names") %>% full_join(shapiro_test, by =
                                                                                          "column_names")
    
    head(shapiro_test)
    
    ### Formatting of the merged file
    all_merged_better = all_merged %>% select(
      column_names,
      useful_features,
      is_numeric,
      useful_features_np,
      feature_completeness.....t..,
      p_values_levene_test,
      p_values_shapiro_test
    ) %>% rename(missingness_percentage = feature_completeness.....t..,
                 feature_names = column_names)
    
    all_merged_better$p_values_l_test_intervals = cut(all_merged_better$p_values_levene_test,
                                                      c(0, 0.001, 0.01, 0.04, 0.06, 1),
                                                      right = FALSE)
    
    all_merged_better$p_values_s_test_intervals = cut(all_merged_better$p_values_shapiro_test,
                                                      c(0, 0.001, 0.01, 0.04, 0.06, 1),
                                                      right = FALSE)
    
    levels(all_merged_better$p_values_l_test_intervals) <-
      c(
        "Significant variation ***",
        "Significant variation **",
        "Significant variation *",
        "Neutral",
        "No significant variation"
      )
    
    levels(all_merged_better$p_values_s_test_intervals) <-
      c(
        "Significantly not normal ***",
        "Significantly not normal **",
        "Significantly not normal *",
        "Neutral",
        "Normal"
      )
    
    levels(all_merged_better$p_values_s_test_intervals)
    
    table(all_merged_better$is_numeric)
    
    #Added in if clause as what if no feature has a missing value or NA
    if(sum(is.na(all_merged_better$is_numeric)) > 0){
      all_merged_better[is.na(all_merged_better$is_numeric), ]$is_numeric <-
      0
    }
    
    all_merged_better$p_values_levene_test = NULL
    all_merged_better$p_values_shapiro_test = NULL
    
    all_merged_better <- 
    
    return(all_merged_better)
  }

# Define UI for application that draws a histogram
ui <- fluidPage(
  # Application title
  titlePanel("Univariate Analysis for Numerical Features"),
  
  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      sliderInput(
        "miss_perc_all",
        "Missingness Percentage Allowed:",
        min = 1,
        max = 100,
        value = 50
      ),
      numericInput("samp_size", "Sample Size for each target level:", 2500),
      numericInput("iterations", "Number of Iterations for tests: ", 1000),
      numericInput("tar_pos", "Target positive value: ", "1"),
      numericInput("tar_neg", "Target negative value: ", "0"),
      numericInput(
        "pval_sig",
        "P-value to determinine significance: ",
        min = 0,
        max = 1,
        step = 0.01,
        0.04
      ),
      numericInput(
        "pval_nsig",
        "P-value cut off for accepting null hypothesis: ",
        min = 0,
        max = 1,
        step = 0.01,
        0.06
      ),
      
      fileInput("file1", "Input a CSV file only: ", accept = c(".csv")),
      actionButton("goButton", "GO!"),
      downloadButton("downloadData", "Download Output"),
      width = 2
    ),
    
    
    # Show a output of the feature analysis
    mainPanel(tableOutput("f"))
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  #Constants
  
  MISSINGNESS_PERCENTAGE = eventReactive(input$goButton, {
    input$miss_perc_all
  })
  ITERATIONS = eventReactive(input$goButton, {
    input$iterations
  })
  TARGET_POSITIVE = eventReactive(input$goButton, {
    input$tar_pos
  })
  TARGET_NEGATIVE = eventReactive(input$goButton, {
    input$tar_neg
  })
  INPUT_FEATURE_FILE = eventReactive(input$goButton, {
    input$file1
  }) #"Features on outreach data MK 2019-06-06.csv"
  SAMPLE_SIZE = eventReactive(input$goButton, {
    input$samp_size
  })
  OUTPUT_FEATURE_FILE = "Formatted Feature Analysis of numeric variables on outreach data.csv"
  P_VALUE_LOWER = eventReactive(input$goButton, {
    input$pval_sig
  })
  P_VALUE_UPPER = eventReactive(input$goButton, {
    input$pval_nsig
  })
  
  output$test <- renderText({
    "Beginning to run"
  })
  #output$f<-renderTable({INPUT_FEATURE_FILE()})
  output$test2 <- renderText({
    "Ended"
  })
  
  import_file <- eventReactive(input$goButton, {
    if (is.null(INPUT_FEATURE_FILE()))
      return(NULL)
    read.csv(INPUT_FEATURE_FILE()$datapath)
  })
  
  datatable <- reactive({
    df <- import_file()
    feature_analysis <-
      fun(
        df,
        MISSINGNESS_PERCENTAGE(),
        ITERATIONS(),
        TARGET_POSITIVE(),
        TARGET_NEGATIVE(),
        SAMPLE_SIZE(),
        P_VALUE_LOWER(),
        P_VALUE_UPPER()
      )
    return(feature_analysis)
  })
  
  output$f <- renderTable({
    datatable()
  })
  
  # Downloadable csv of selected dataset ----
  output$downloadData <- downloadHandler(
    filename = "feature_analysis.csv",
    content = function(file) {
      write.csv(datatable(), file, row.names = FALSE)
    }
  )
  
  
}

# Run the application
shinyApp(ui = ui, server = server)
  