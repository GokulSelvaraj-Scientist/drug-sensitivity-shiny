# ============================================================
# Drug Sensitivity Prediction: Interactive Shiny App
# Author: Gokul Selvaraj
# GitHub: GokulSelvaraj-Scientist
# Description: Interactive dashboard for exploring and predicting
#              drug sensitivity in cancer cell lines using
#              simulated GDSC-style data
# ============================================================

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(tidyr)
library(randomForest)
library(caret)
library(pROC)
library(RColorBrewer)
library(DT)

# ============================================================
# SIMULATE REALISTIC GDSC-STYLE DATASET
# Based on published distributions from the Genomics of Drug
# Sensitivity in Cancer (GDSC) database
# ============================================================

set.seed(42)
n_lines <- 200

# Cancer types
cancer_types <- c("Breast", "Lung", "Colorectal", "Leukemia", "Melanoma",
                  "Ovarian", "Pancreatic", "Prostate")

# Drugs with known mechanisms
drugs <- c("Erlotinib", "Lapatinib", "Sorafenib", "Imatinib",
           "Docetaxel", "Gemcitabine", "Cisplatin", "Vemurafenib")

# Simulate genomic features
cell_lines <- data.frame(
  CellLine     = paste0("CL_", seq_len(n_lines)),
  CancerType   = sample(cancer_types, n_lines, replace = TRUE),
  TP53_mut     = rbinom(n_lines, 1, 0.4),
  KRAS_mut     = rbinom(n_lines, 1, 0.25),
  BRAF_mut     = rbinom(n_lines, 1, 0.15),
  EGFR_mut     = rbinom(n_lines, 1, 0.2),
  PIK3CA_mut   = rbinom(n_lines, 1, 0.3),
  MYC_amp      = rbinom(n_lines, 1, 0.2),
  PTEN_loss    = rbinom(n_lines, 1, 0.25),
  MSI_status   = rbinom(n_lines, 1, 0.15),
  TMB          = rnbinom(n_lines, mu = 5, size = 2),
  Ploidy       = round(rnorm(n_lines, mean = 2.5, sd = 0.8), 1),
  DoubleTime   = round(rnorm(n_lines, mean = 30, sd = 10), 1)
)

# Simulate IC50 values with biologically realistic drug-mutation associations
simulate_ic50 <- function(cell_lines, drug) {
  base_ic50 <- rnorm(nrow(cell_lines), mean = 2, sd = 1.5)

  ic50 <- switch(drug,
    "Erlotinib"   = base_ic50 - 2 * cell_lines$EGFR_mut + 1.5 * cell_lines$KRAS_mut,
    "Lapatinib"   = base_ic50 - 1.8 * cell_lines$EGFR_mut + 1.2 * cell_lines$KRAS_mut,
    "Sorafenib"   = base_ic50 - 1.5 * cell_lines$BRAF_mut + rnorm(nrow(cell_lines), 0, 0.5),
    "Imatinib"    = base_ic50 - 2.5 * (cell_lines$CancerType == "Leukemia"),
    "Docetaxel"   = base_ic50 + 0.5 * cell_lines$TP53_mut + rnorm(nrow(cell_lines), 0, 0.8),
    "Gemcitabine" = base_ic50 - 0.8 * cell_lines$KRAS_mut + rnorm(nrow(cell_lines), 0, 0.6),
    "Cisplatin"   = base_ic50 - 1.2 * cell_lines$PTEN_loss + rnorm(nrow(cell_lines), 0, 0.7),
    "Vemurafenib" = base_ic50 - 3.0 * cell_lines$BRAF_mut + 2.0 * cell_lines$KRAS_mut,
    base_ic50
  )
  ic50 + rnorm(nrow(cell_lines), 0, 0.3)
}

# Create full sensitivity matrix
sensitivity_data <- lapply(drugs, function(drug) {
  data.frame(
    CellLine  = cell_lines$CellLine,
    Drug      = drug,
    IC50_log  = simulate_ic50(cell_lines, drug),
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

# Merge with cell line features
full_data <- sensitivity_data %>%
  left_join(cell_lines, by = "CellLine") %>%
  mutate(
    Sensitive = factor(ifelse(IC50_log < median(IC50_log), "Sensitive", "Resistant"),
                       levels = c("Sensitive", "Resistant"))
  )

# ============================================================
# UI
# ============================================================

ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(
    title = "Drug Sensitivity Explorer",
    titleWidth = 280
  ),

  dashboardSidebar(
    width = 280,
    sidebarMenu(
      menuItem("Overview",         tabName = "overview",    icon = icon("chart-bar")),
      menuItem("Drug Explorer",    tabName = "explorer",    icon = icon("flask")),
      menuItem("ML Prediction",    tabName = "prediction",  icon = icon("brain")),
      menuItem("Biomarker Analysis", tabName = "biomarker", icon = icon("dna"))
    ),
    hr(),
    # Global filters
    selectInput("drug_select", "Select Drug:",
                choices  = drugs,
                selected = "Erlotinib"),
    selectInput("cancer_select", "Cancer Type:",
                choices  = c("All", cancer_types),
                selected = "All"),
    hr(),
    p(style = "padding-left:15px; color:#aaa; font-size:11px;",
      "Data: Simulated GDSC-style dataset\nn=200 cancer cell lines\n8 drugs | 12 genomic features")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #f4f6f9; }
      .box { border-radius: 8px; }
      .info-box { border-radius: 8px; }
    "))),

    tabItems(

      # ---- TAB 1: OVERVIEW ----
      tabItem(tabName = "overview",
        fluidRow(
          infoBoxOutput("n_celllines", width = 3),
          infoBoxOutput("n_drugs",     width = 3),
          infoBoxOutput("n_features",  width = 3),
          infoBoxOutput("n_cancers",   width = 3)
        ),
        fluidRow(
          box(title = "IC50 Distribution by Cancer Type", width = 6, status = "primary",
              plotOutput("overview_boxplot", height = "350px")),
          box(title = "Sensitivity Heatmap: All Drugs", width = 6, status = "primary",
              plotOutput("sensitivity_heatmap", height = "350px"))
        ),
        fluidRow(
          box(title = "Mutation Frequency by Cancer Type", width = 12, status = "info",
              plotOutput("mutation_freq", height = "300px"))
        )
      ),

      # ---- TAB 2: DRUG EXPLORER ----
      tabItem(tabName = "explorer",
        fluidRow(
          box(title = "IC50 Distribution", width = 6, status = "primary",
              plotOutput("ic50_dist", height = "300px")),
          box(title = "IC50 by Mutation Status", width = 6, status = "primary",
              selectInput("mut_compare", "Compare by mutation:",
                          choices = c("TP53_mut", "KRAS_mut", "BRAF_mut",
                                      "EGFR_mut", "PIK3CA_mut", "PTEN_loss")),
              plotOutput("ic50_by_mutation", height = "260px"))
        ),
        fluidRow(
          box(title = "Sensitivity vs Genomic Features", width = 6, status = "info",
              selectInput("scatter_x", "X axis feature:",
                          choices = c("TMB", "Ploidy", "DoubleTime")),
              plotOutput("scatter_plot", height = "280px")),
          box(title = "Drug Sensitivity Summary Table", width = 6, status = "info",
              DTOutput("summary_table"))
        )
      ),

      # ---- TAB 3: ML PREDICTION ----
      tabItem(tabName = "prediction",
        fluidRow(
          box(title = "Model Training Controls", width = 4, status = "warning",
              sliderInput("train_split", "Training set %:", min = 60, max = 90,
                          value = 80, step = 5),
              selectInput("ml_drug", "Drug to predict:",
                          choices = drugs, selected = "Erlotinib"),
              selectInput("ml_model", "Model:",
                          choices = c("Random Forest" = "rf",
                                      "Logistic Regression" = "glm")),
              actionButton("train_model", "Train Model",
                           class = "btn-primary btn-block",
                           icon  = icon("play")),
              hr(),
              verbatimTextOutput("model_summary")
          ),
          box(title = "ROC Curve", width = 4, status = "primary",
              plotOutput("roc_plot", height = "320px")),
          box(title = "Confusion Matrix", width = 4, status = "primary",
              plotOutput("cm_plot", height = "320px"))
        ),
        fluidRow(
          box(title = "Feature Importance", width = 12, status = "info",
              plotOutput("importance_plot", height = "280px"))
        )
      ),

      # ---- TAB 4: BIOMARKER ANALYSIS ----
      tabItem(tabName = "biomarker",
        fluidRow(
          box(title = "Biomarker Selection", width = 4, status = "warning",
              selectInput("biomarker_drug", "Drug:", choices = drugs),
              selectInput("biomarker_mut", "Mutation/Feature:",
                          choices = c("TP53_mut", "KRAS_mut", "BRAF_mut",
                                      "EGFR_mut", "PIK3CA_mut", "PTEN_loss",
                                      "MSI_status", "MYC_amp")),
              hr(),
              p("Statistical test: Wilcoxon rank-sum test"),
              verbatimTextOutput("stat_test_result")
          ),
          box(title = "IC50 by Biomarker Status", width = 8, status = "primary",
              plotOutput("biomarker_plot", height = "380px"))
        ),
        fluidRow(
          box(title = "Biomarker Summary Across All Drugs", width = 12, status = "info",
              plotOutput("biomarker_heatmap", height = "320px"))
        )
      )
    )
  )
)

# ============================================================
# SERVER
# ============================================================

server <- function(input, output, session) {

  # --- Reactive filtered data ---
  filtered_data <- reactive({
    df <- full_data %>% filter(Drug == input$drug_select)
    if (input$cancer_select != "All") {
      df <- df %>% filter(CancerType == input$cancer_select)
    }
    df
  })

  # ---- OVERVIEW TAB ----
  output$n_celllines <- renderInfoBox({
    infoBox("Cell Lines", 200, icon = icon("vial"), color = "blue")
  })
  output$n_drugs <- renderInfoBox({
    infoBox("Drugs", 8, icon = icon("pills"), color = "green")
  })
  output$n_features <- renderInfoBox({
    infoBox("Genomic Features", 12, icon = icon("dna"), color = "purple")
  })
  output$n_cancers <- renderInfoBox({
    infoBox("Cancer Types", 8, icon = icon("hospital"), color = "red")
  })

  output$overview_boxplot <- renderPlot({
    full_data %>%
      filter(Drug == input$drug_select) %>%
      ggplot(aes(x = reorder(CancerType, IC50_log, median), y = IC50_log, fill = CancerType)) +
      geom_boxplot(alpha = 0.7, outlier.size = 1) +
      coord_flip() +
      scale_fill_brewer(palette = "Set2") +
      labs(x = "Cancer Type", y = "Log IC50", title = paste("IC50 Distribution —", input$drug_select)) +
      theme_classic(base_size = 12) +
      theme(legend.position = "none", plot.title = element_text(face = "bold"))
  })

  output$sensitivity_heatmap <- renderPlot({
    heatmap_data <- full_data %>%
      group_by(Drug, CancerType) %>%
      summarise(Mean_IC50 = mean(IC50_log), .groups = "drop")

    ggplot(heatmap_data, aes(x = CancerType, y = Drug, fill = Mean_IC50)) +
      geom_tile(color = "white") +
      scale_fill_gradient2(low = "#2A9D8F", mid = "white", high = "#E76F51",
                           midpoint = median(heatmap_data$Mean_IC50)) +
      labs(x = "Cancer Type", y = "Drug", fill = "Mean\nLog IC50",
           title = "Drug Sensitivity Landscape") +
      theme_classic(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(face = "bold"))
  })

  output$mutation_freq <- renderPlot({
    mut_data <- cell_lines %>%
      group_by(CancerType) %>%
      summarise(across(c(TP53_mut, KRAS_mut, BRAF_mut, EGFR_mut, PIK3CA_mut, PTEN_loss),
                       mean, .names = "{.col}")) %>%
      pivot_longer(-CancerType, names_to = "Mutation", values_to = "Frequency")

    ggplot(mut_data, aes(x = CancerType, y = Frequency, fill = Mutation)) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_fill_brewer(palette = "Set1") +
      labs(x = "Cancer Type", y = "Mutation Frequency", title = "Mutation Frequencies by Cancer Type") +
      theme_classic(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(face = "bold"))
  })

  # ---- DRUG EXPLORER TAB ----
  output$ic50_dist <- renderPlot({
    ggplot(filtered_data(), aes(x = IC50_log, fill = Sensitive)) +
      geom_histogram(bins = 30, alpha = 0.7, color = "white") +
      scale_fill_manual(values = c("Sensitive" = "#2A9D8F", "Resistant" = "#E76F51")) +
      geom_vline(xintercept = median(filtered_data()$IC50_log),
                 linetype = "dashed", color = "black") +
      labs(x = "Log IC50", y = "Count", fill = "Status",
           title = paste("IC50 Distribution —", input$drug_select),
           subtitle = "Dashed line = median threshold") +
      theme_classic(base_size = 12) +
      theme(plot.title = element_text(face = "bold"))
  })

  output$ic50_by_mutation <- renderPlot({
    df <- filtered_data() %>%
      mutate(MutStatus = factor(ifelse(.data[[input$mut_compare]] == 1,
                                       "Mutant", "Wild-type")))

    ggplot(df, aes(x = MutStatus, y = IC50_log, fill = MutStatus)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      geom_jitter(width = 0.2, size = 1.5, alpha = 0.5) +
      scale_fill_manual(values = c("Mutant" = "#E76F51", "Wild-type" = "#457B9D")) +
      labs(x = input$mut_compare, y = "Log IC50",
           title = paste(input$drug_select, "sensitivity by", input$mut_compare)) +
      theme_classic(base_size = 12) +
      theme(legend.position = "none", plot.title = element_text(face = "bold"))
  })

  output$scatter_plot <- renderPlot({
    ggplot(filtered_data(), aes_string(x = input$scatter_x, y = "IC50_log", color = "Sensitive")) +
      geom_point(size = 2, alpha = 0.6) +
      geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
      scale_color_manual(values = c("Sensitive" = "#2A9D8F", "Resistant" = "#E76F51")) +
      labs(y = "Log IC50", color = "Status",
           title = paste("IC50 vs", input$scatter_x)) +
      theme_classic(base_size = 12) +
      theme(plot.title = element_text(face = "bold"))
  })

  output$summary_table <- renderDT({
    filtered_data() %>%
      group_by(CancerType) %>%
      summarise(
        N         = n(),
        Mean_IC50 = round(mean(IC50_log), 2),
        SD_IC50   = round(sd(IC50_log), 2),
        Pct_Sensitive = round(mean(Sensitive == "Sensitive") * 100, 1)
      ) %>%
      arrange(Mean_IC50) %>%
      datatable(options = list(pageLength = 8, dom = "t"),
                colnames = c("Cancer Type", "N", "Mean Log IC50", "SD", "% Sensitive"))
  })

  # ---- ML PREDICTION TAB ----
  model_results <- eventReactive(input$train_model, {
    withProgress(message = "Training model...", {
      df <- full_data %>%
        filter(Drug == input$ml_drug) %>%
        select(Sensitive, TP53_mut, KRAS_mut, BRAF_mut, EGFR_mut,
               PIK3CA_mut, MYC_amp, PTEN_loss, MSI_status, TMB, Ploidy, DoubleTime) %>%
        mutate(across(c(TP53_mut, KRAS_mut, BRAF_mut, EGFR_mut,
                        PIK3CA_mut, MYC_amp, PTEN_loss, MSI_status), factor))

      set.seed(42)
      train_idx  <- createDataPartition(df$Sensitive, p = input$train_split / 100, list = FALSE)
      train_data <- df[train_idx, ]
      test_data  <- df[-train_idx, ]

      ctrl <- trainControl(method = "cv", number = 5,
                           classProbs = TRUE, summaryFunction = twoClassSummary)

      if (input$ml_model == "rf") {
        model <- train(Sensitive ~ ., data = train_data, method = "rf",
                       trControl = ctrl, metric = "ROC")
      } else {
        model <- train(Sensitive ~ ., data = train_data, method = "glm",
                       family = "binomial", trControl = ctrl, metric = "ROC")
      }

      preds  <- predict(model, test_data)
      probs  <- predict(model, test_data, type = "prob")
      cm     <- confusionMatrix(preds, test_data$Sensitive, positive = "Sensitive")
      roc_obj <- roc(test_data$Sensitive, probs$Sensitive,
                     levels = c("Resistant", "Sensitive"))

      list(model = model, cm = cm, roc = roc_obj, test_data = test_data, probs = probs)
    })
  })

  output$model_summary <- renderPrint({
    req(model_results())
    cm <- model_results()$cm
    cat("Accuracy:   ", round(cm$overall["Accuracy"], 3), "\n")
    cat("Sensitivity:", round(cm$byClass["Sensitivity"], 3), "\n")
    cat("Specificity:", round(cm$byClass["Specificity"], 3), "\n")
    cat("AUC:        ", round(auc(model_results()$roc), 3), "\n")
  })

  output$roc_plot <- renderPlot({
    req(model_results())
    roc_obj <- model_results()$roc
    roc_df  <- data.frame(FPR = 1 - roc_obj$specificities, TPR = roc_obj$sensitivities)

    ggplot(roc_df, aes(x = FPR, y = TPR)) +
      geom_line(color = "#457B9D", linewidth = 1.5) +
      geom_abline(linetype = "dashed", color = "grey50") +
      annotate("text", x = 0.6, y = 0.2,
               label = paste0("AUC = ", round(auc(roc_obj), 3)),
               size = 5, color = "#457B9D", fontface = "bold") +
      labs(title = "ROC Curve", x = "False Positive Rate", y = "True Positive Rate") +
      theme_classic(base_size = 12) +
      theme(plot.title = element_text(face = "bold"))
  })

  output$cm_plot <- renderPlot({
    req(model_results())
    cm_df <- as.data.frame(model_results()$cm$table) %>%
      rename(Predicted = Prediction, Actual = Reference)

    ggplot(cm_df, aes(x = Actual, y = Predicted, fill = Freq)) +
      geom_tile(color = "white") +
      geom_text(aes(label = Freq), size = 8, fontface = "bold") +
      scale_fill_gradient(low = "#F1FAEE", high = "#E63946") +
      labs(title = "Confusion Matrix",
           subtitle = paste0("Accuracy: ",
                             round(model_results()$cm$overall["Accuracy"] * 100, 1), "%")) +
      theme_classic(base_size = 12) +
      theme(plot.title = element_text(face = "bold"))
  })

  output$importance_plot <- renderPlot({
    req(model_results())
    imp <- varImp(model_results()$model)$importance %>%
      tibble::rownames_to_column("Feature") %>%
      arrange(desc(Overall)) %>%
      head(10)

    ggplot(imp, aes(x = reorder(Feature, Overall), y = Overall, fill = Overall)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      scale_fill_gradient(low = "#A8DADC", high = "#E63946") +
      labs(title = "Feature Importance", x = "Feature", y = "Importance Score") +
      theme_classic(base_size = 12) +
      theme(plot.title = element_text(face = "bold"), legend.position = "none")
  })

  # ---- BIOMARKER TAB ----
  output$biomarker_plot <- renderPlot({
    df <- full_data %>%
      filter(Drug == input$biomarker_drug) %>%
      mutate(BiomarkerStatus = factor(ifelse(.data[[input$biomarker_mut]] == 1,
                                             "Mutant", "Wild-type")))

    ggplot(df, aes(x = BiomarkerStatus, y = IC50_log, fill = BiomarkerStatus)) +
      geom_violin(alpha = 0.6, trim = FALSE) +
      geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
      geom_jitter(width = 0.1, size = 1.5, alpha = 0.4) +
      scale_fill_manual(values = c("Mutant" = "#E76F51", "Wild-type" = "#2A9D8F")) +
      facet_wrap(~ CancerType, ncol = 4) +
      labs(x = input$biomarker_mut, y = "Log IC50",
           title = paste(input$biomarker_drug, "sensitivity by", input$biomarker_mut),
           subtitle = "Stratified by cancer type") +
      theme_classic(base_size = 11) +
      theme(legend.position = "none", plot.title = element_text(face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1))
  })

  output$stat_test_result <- renderPrint({
    df <- full_data %>%
      filter(Drug == input$biomarker_drug)

    mutant   <- df$IC50_log[df[[input$biomarker_mut]] == 1]
    wildtype <- df$IC50_log[df[[input$biomarker_mut]] == 0]

    if (length(mutant) < 3 || length(wildtype) < 3) {
      cat("Insufficient data for test\n")
      return()
    }

    test <- wilcox.test(mutant, wildtype)
    cat("Wilcoxon rank-sum test\n")
    cat("Mutant n =", length(mutant), "| Wild-type n =", length(wildtype), "\n")
    cat("Mutant median IC50:", round(median(mutant), 3), "\n")
    cat("Wild-type median IC50:", round(median(wildtype), 3), "\n")
    cat("W =", test$statistic, "\n")
    cat("p-value =", round(test$p.value, 4), "\n")
    cat(ifelse(test$p.value < 0.05, "*** Significant difference", "Not significant"), "\n")
  })

  output$biomarker_heatmap <- renderPlot({
    mutations <- c("TP53_mut", "KRAS_mut", "BRAF_mut", "EGFR_mut",
                   "PIK3CA_mut", "PTEN_loss", "MSI_status", "MYC_amp")

    effect_data <- expand.grid(Drug = drugs, Mutation = mutations) %>%
      rowwise() %>%
      mutate(
        Effect = {
          df    <- full_data %>% filter(Drug == Drug)
          mut   <- df$IC50_log[df[[as.character(Mutation)]] == 1]
          wt    <- df$IC50_log[df[[as.character(Mutation)]] == 0]
          if (length(mut) < 3 || length(wt) < 3) NA
          else median(mut) - median(wt)
        }
      ) %>%
      ungroup()

    ggplot(effect_data, aes(x = Mutation, y = Drug, fill = Effect)) +
      geom_tile(color = "white") +
      scale_fill_gradient2(low = "#2A9D8F", mid = "white", high = "#E76F51",
                           midpoint = 0, na.value = "grey80") +
      labs(x = "Mutation", y = "Drug", fill = "Median IC50\nDifference\n(Mut - WT)",
           title = "Biomarker Effect Map: Mutation Impact on Drug Sensitivity") +
      theme_classic(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title  = element_text(face = "bold"))
  })
}

# ============================================================
# RUN APP
# ============================================================
shinyApp(ui = ui, server = server)
